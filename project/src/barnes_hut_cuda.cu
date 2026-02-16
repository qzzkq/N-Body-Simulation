#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include "barnes_hut_cuda.cuh"
#include "object.hpp"

#define THREADS_PER_BLOCK 256
#define MAX_STACK 64
#define THETA 0.5
#define G_CONST 6.6743e-20


#define SOFTENING_SQ 0.1 

struct __align__(16) GpuNode {
    double3 center;
    double mass;
    double3 bMin;
    double size;
    int children[8];
    int bodyIdx;
};

static double4* d_posMass = nullptr;
static double3* d_vel = nullptr;
static double3* d_acc = nullptr;
static GpuNode* d_nodes = nullptr;
static int* d_rootIdx = nullptr;

static std::vector<GpuNode> hostNodes;
static size_t currentCapacity = 0;

__global__ void computeForcesKernel(
    const double4* __restrict__ posMass,
    double3* __restrict__ vel,
    double3* __restrict__ acc,
    const GpuNode* __restrict__ nodes,
    int numBodies,
    int rootIdx,
    double dt
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numBodies) return;

    double4 myP = posMass[idx];
    double3 myPos = make_double3(myP.x, myP.y, myP.z);
    double3 newForce = make_double3(0.0, 0.0, 0.0);

    int stack[MAX_STACK];
    int stackTop = 0;
    
    stack[stackTop++] = rootIdx;

    while (stackTop > 0) {
        int nodeIdx = stack[--stackTop];
        GpuNode n = nodes[nodeIdx];

        double dx = n.center.x - myPos.x;
        double dy = n.center.y - myPos.y;
        double dz = n.center.z - myPos.z;
        
        double distSq = dx*dx + dy*dy + dz*dz + SOFTENING_SQ; 
        double dist = sqrt(distSq);

        bool isLeaf = (n.bodyIdx != -1);
        
        if (isLeaf && n.bodyIdx == idx) continue;

        if (isLeaf || (n.size / dist < THETA)) {
            double f = (G_CONST * n.mass) / (distSq * dist);
            newForce.x += f * dx;
            newForce.y += f * dy;
            newForce.z += f * dz;
        } else {
            #pragma unroll
            for (int i = 0; i < 8; ++i) {
                if (n.children[i] != -1) {
                    if (stackTop < MAX_STACK) {
                        stack[stackTop++] = n.children[i];
                    }
                }
            }
        }
    }

    vel[idx].x += newForce.x * dt;
    vel[idx].y += newForce.y * dt;
    vel[idx].z += newForce.z * dt;
}

__global__ void integratePositionKernel(
    double4* posMass,
    const double3* vel,
    int numBodies,
    double dt
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numBodies) return;

    posMass[idx].x += vel[idx].x * dt;
    posMass[idx].y += vel[idx].y * dt;
    posMass[idx].z += vel[idx].z * dt;
}

struct BuildNode {
    double mass;
    double cx, cy, cz;
    double x, y, z, size;
    int children[8];
    int bodyIdx;
    
    void init(double _x, double _y, double _z, double _s) {
        mass = 0; cx=0; cy=0; cz=0;
        x=_x; y=_y; z=_z; size=_s;
        bodyIdx = -1;
        for(int i=0; i<8; ++i) children[i] = -1;
    }
};

static std::vector<BuildNode> buildPool;
static int nodesUsed = 0;

static int newBuildNode(double x, double y, double z, double s) {
    if (nodesUsed >= buildPool.size()) {
        buildPool.resize(nodesUsed + 4096);
    }
    int idx = nodesUsed++;
    buildPool[idx].init(x, y, z, s);
    return idx;
}

static int getOctant(double tx, double ty, double tz, double size, double px, double py, double pz) {
    int idx = 0;
    if (px >= tx) idx |= 1;
    if (py >= ty) idx |= 2;
    if (pz >= tz) idx |= 4;
    return idx;
}

static void insertBody(int nodeIdx, int bodyIdx, const glm::dvec3& pos, double mass) {
    BuildNode* node = &buildPool[nodeIdx];

    if (node->size < 1e-9) { 
        return;
    }

    bool hasChildren = false;
    for(int i = 0; i < 8; ++i) {
        if(node->children[i] != -1) {
            hasChildren = true;
            break;
        }
    }

    if (node->mass == 0.0 && node->bodyIdx == -1 && !hasChildren) {
        node->bodyIdx = bodyIdx;
        node->mass = mass;
        node->cx = pos.x; node->cy = pos.y; node->cz = pos.z;
        return;
    }

    if (node->bodyIdx != -1) {
        int oldBody = node->bodyIdx;
        node->bodyIdx = -1; 
        
        int oct = getOctant(node->x, node->y, node->z, node->size, node->cx, node->cy, node->cz);
        
        if (node->children[oct] == -1) {
            double qs = node->size * 0.25;
            double nx = node->x + ((oct&1) ? qs : -qs);
            double ny = node->y + ((oct&2) ? qs : -qs);
            double nz = node->z + ((oct&4) ? qs : -qs);
            
            int childIdx = newBuildNode(nx, ny, nz, node->size * 0.5);
            node = &buildPool[nodeIdx]; 
            node->children[oct] = childIdx;
        }

        insertBody(node->children[oct], oldBody, glm::dvec3(node->cx, node->cy, node->cz), node->mass);
        
        node = &buildPool[nodeIdx]; 

        node->mass = 0;
        node->cx = 0; node->cy = 0; node->cz = 0;
    }

    int oct = getOctant(node->x, node->y, node->z, node->size, pos.x, pos.y, pos.z);
    
    if (node->children[oct] == -1) {
        double qs = node->size * 0.25;
        double nx = node->x + ((oct&1) ? qs : -qs);
        double ny = node->y + ((oct&2) ? qs : -qs);
        double nz = node->z + ((oct&4) ? qs : -qs);
        
        int childIdx = newBuildNode(nx, ny, nz, node->size * 0.5);
        node = &buildPool[nodeIdx]; 
        node->children[oct] = childIdx;
    }
    
    insertBody(node->children[oct], bodyIdx, pos, mass);
}

static void computeMassDist(int nodeIdx) {
    BuildNode* node = &buildPool[nodeIdx];
    if (node->bodyIdx != -1) return;

    double m = 0, cx = 0, cy = 0, cz = 0;
    for (int i = 0; i < 8; ++i) {
        if (node->children[i] != -1) {
            computeMassDist(node->children[i]);
            BuildNode& child = buildPool[node->children[i]];
            m += child.mass;
            cx += child.cx * child.mass;
            cy += child.cy * child.mass;
            cz += child.cz * child.mass;
        }
    }
    node->mass = m;
    if (m > 0) {
        node->cx = cx / m;
        node->cy = cy / m;
        node->cz = cz / m;
    }
}

static int flattenTree(int nodeIdx, int& outIdx) {
    int currentIdx = outIdx++;
    BuildNode& bn = buildPool[nodeIdx];
    GpuNode& gn = hostNodes[currentIdx];

    gn.center = make_double3(bn.cx, bn.cy, bn.cz);
    gn.mass = bn.mass;
    gn.bMin = make_double3(bn.x, bn.y, bn.z);
    gn.size = bn.size;
    gn.bodyIdx = bn.bodyIdx;

    for (int i = 0; i < 8; ++i) {
        if (bn.children[i] != -1) {
            gn.children[i] = flattenTree(bn.children[i], outIdx);
        } else {
            gn.children[i] = -1;
        }
    }
    return currentIdx;
}

void InitBarnesHutCUDA(size_t maxObjects) {
    currentCapacity = maxObjects + 4096;
    cudaMalloc(&d_posMass, currentCapacity * sizeof(double4));
    cudaMalloc(&d_vel, currentCapacity * sizeof(double3));
    cudaMalloc(&d_acc, currentCapacity * sizeof(double3));
    cudaMalloc(&d_nodes, currentCapacity * 4 * sizeof(GpuNode));
    cudaMalloc(&d_rootIdx, sizeof(int));
    
    hostNodes.resize(currentCapacity * 4);
    buildPool.reserve(currentCapacity * 4);
}

void CleanupBarnesHutCUDA() {
    if (d_posMass) cudaFree(d_posMass);
    if (d_vel) cudaFree(d_vel);
    if (d_acc) cudaFree(d_acc);
    if (d_nodes) cudaFree(d_nodes);
    if (d_rootIdx) cudaFree(d_rootIdx);
}

void simulationStepBarnesHutCUDA(std::vector<Object>& objs, float dt, bool pause, int iterations) {
    if (objs.empty() || pause) return;
    if (iterations < 1) iterations = 1;

    size_t N = objs.size();
    if (N > currentCapacity) {
        CleanupBarnesHutCUDA();
        InitBarnesHutCUDA(N * 2);
    }

    std::vector<double4> h_posMass(N);
    std::vector<double3> h_vel(N);
    
    for (size_t i = 0; i < N; ++i) {
        auto p = objs[i].GetPos();
        h_posMass[i] = make_double4(p.x, p.y, p.z, objs[i].mass);
        h_vel[i] = make_double3(objs[i].velocity.x, objs[i].velocity.y, objs[i].velocity.z);
    }

    cudaMemcpy(d_posMass, h_posMass.data(), N * sizeof(double4), cudaMemcpyHostToDevice);
    cudaMemcpy(d_vel, h_vel.data(), N * sizeof(double3), cudaMemcpyHostToDevice);

    for (int iter = 0; iter < iterations; ++iter) {
        if (iter > 0) {
            cudaMemcpy(h_posMass.data(), d_posMass, N * sizeof(double4), cudaMemcpyDeviceToHost);
        }

        double minX = 1e30, maxX = -1e30;
        double minY = 1e30, maxY = -1e30;
        double minZ = 1e30, maxZ = -1e30;

        for (size_t i = 0; i < N; ++i) {
            double x = h_posMass[i].x;
            double y = h_posMass[i].y;
            double z = h_posMass[i].z;

            if (x < minX) minX = x; if (x > maxX) maxX = x;
            if (y < minY) minY = y; if (y > maxY) maxY = y;
            if (z < minZ) minZ = z; if (z > maxZ) maxZ = z;
        }

        nodesUsed = 0;
        double maxSpan = std::max({maxX - minX, maxY - minY, maxZ - minZ});
        maxSpan *= 1.01; 
        
        double cx = (minX + maxX) * 0.5;
        double cy = (minY + maxY) * 0.5;
        double cz = (minZ + maxZ) * 0.5;

        int rootBuildIdx = newBuildNode(cx, cy, cz, maxSpan);

        for (int i = 0; i < N; ++i) {
            glm::dvec3 pos(h_posMass[i].x, h_posMass[i].y, h_posMass[i].z);
            double mass = h_posMass[i].w;
            insertBody(rootBuildIdx, i, pos, mass);
        }
        computeMassDist(rootBuildIdx);

        int outIdx = 0;
        int rootGpuIdx = flattenTree(rootBuildIdx, outIdx);

        cudaMemcpy(d_nodes, hostNodes.data(), outIdx * sizeof(GpuNode), cudaMemcpyHostToDevice);

        int blocks = (N + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
        computeForcesKernel<<<blocks, THREADS_PER_BLOCK>>>(
            d_posMass, d_vel, d_acc, d_nodes, N, rootGpuIdx, (double)dt
        );
        
        integratePositionKernel<<<blocks, THREADS_PER_BLOCK>>>(
            d_posMass, d_vel, N, (double)dt
        );
        
    }

    cudaMemcpy(h_posMass.data(), d_posMass, N * sizeof(double4), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_vel.data(), d_vel, N * sizeof(double3), cudaMemcpyDeviceToHost);

    for (size_t i = 0; i < N; ++i) {
        objs[i].position.x = h_posMass[i].x;
        objs[i].position.y = h_posMass[i].y;
        objs[i].position.z = h_posMass[i].z;
        objs[i].velocity.x = h_vel[i].x;
        objs[i].velocity.y = h_vel[i].y;
        objs[i].velocity.z = h_vel[i].z;
    }
}