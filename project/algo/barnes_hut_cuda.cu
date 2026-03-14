#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include "barnes_hut_cuda.cuh"
#include "object.hpp"
#include "physics.hpp"

#define THREADS_PER_BLOCK 256
#define MAX_STACK 64
#define THETA 0.5

static constexpr double G_CONST = physics::G;
static constexpr double SOFTENING_AU = 1.0e6 * physics::METERS_TO_AU; // 1000 km
static constexpr double SOFTENING_SQ = SOFTENING_AU * SOFTENING_AU;

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
static int* d_nodeCounter = nullptr; 
static double* d_minMax = nullptr; 

static size_t currentCapacity = 0;

__device__ void atomicMinDouble(double* addr, double value) {
    double old = *addr, assumed;
    if (old <= value) return;
    do {
        assumed = old;
        old = __longlong_as_double(atomicCAS((unsigned long long int*)addr, __double_as_longlong(assumed), __double_as_longlong(value)));
    } while (assumed != old && old > value);
}

__device__ void atomicMaxDouble(double* addr, double value) {
    double old = *addr, assumed;
    if (old >= value) return;
    do {
        assumed = old;
        old = __longlong_as_double(atomicCAS((unsigned long long int*)addr, __double_as_longlong(assumed), __double_as_longlong(value)));
    } while (assumed != old && old < value);
}

__device__ int getOctant(double tx, double ty, double tz, double px, double py, double pz) {
    int idx = 0;
    if (px >= tx) idx |= 1;
    if (py >= ty) idx |= 2;
    if (pz >= tz) idx |= 4;
    return idx;
}

__device__ void initNode(GpuNode* nodes, int idx, double3 bMin, double size) {
    nodes[idx].bMin = bMin;
    nodes[idx].size = size;
    nodes[idx].mass = 0.0;
    nodes[idx].center = make_double3(0.0, 0.0, 0.0);
    nodes[idx].bodyIdx = -1;
    #pragma unroll
    for (int i = 0; i < 8; ++i) {
        nodes[idx].children[i] = -1;
    }
}

__global__ void computeBoundingBoxKernel(const double4* posMass, int numBodies, double* d_minMax) {
    __shared__ double s_minMax[6]; 

    if (threadIdx.x == 0) {
        s_minMax[0] = 1e30; s_minMax[1] = -1e30;
        s_minMax[2] = 1e30; s_minMax[3] = -1e30;
        s_minMax[4] = 1e30; s_minMax[5] = -1e30;
    }
    __syncthreads();

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < numBodies) {
        double4 p = posMass[idx];
        atomicMinDouble(&s_minMax[0], p.x); atomicMaxDouble(&s_minMax[1], p.x);
        atomicMinDouble(&s_minMax[2], p.y); atomicMaxDouble(&s_minMax[3], p.y);
        atomicMinDouble(&s_minMax[4], p.z); atomicMaxDouble(&s_minMax[5], p.z);
    }
    __syncthreads();

    if (threadIdx.x == 0) {
        atomicMinDouble(&d_minMax[0], s_minMax[0]); atomicMaxDouble(&d_minMax[1], s_minMax[1]);
        atomicMinDouble(&d_minMax[2], s_minMax[2]); atomicMaxDouble(&d_minMax[3], s_minMax[3]);
        atomicMinDouble(&d_minMax[4], s_minMax[4]); atomicMaxDouble(&d_minMax[5], s_minMax[5]);
    }
}

__global__ void initRootKernel(GpuNode* nodes, double cx, double cy, double cz, double maxSpan) {
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        initNode(nodes, 0, make_double3(cx, cy, cz), maxSpan);
    }
}

__global__ void buildTreeKernel(const double4* __restrict__ posMass, GpuNode* __restrict__ nodes, int numBodies, int* nodeCounter, int maxNodes) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numBodies) return;

    double4 p = posMass[i];
    double3 pos = make_double3(p.x, p.y, p.z);
    double mass = p.w;
    if (mass <= 0.0) return;

    int curr = 0; 
    
    while (true) {
        atomicAdd(&nodes[curr].mass, mass);
        atomicAdd(&nodes[curr].center.x, pos.x * mass);
        atomicAdd(&nodes[curr].center.y, pos.y * mass);
        atomicAdd(&nodes[curr].center.z, pos.z * mass);
        
        if (nodes[curr].size < 1e-8) return; 
        
        while (true) {
            int oldBody = atomicCAS(&nodes[curr].bodyIdx, -1, i);
            
            if (oldBody == -1 || oldBody == i) return; 
            
            if (oldBody >= 0) {
                if (atomicCAS(&nodes[curr].bodyIdx, oldBody, -2) == oldBody) {
                    double3 currPos = nodes[curr].bMin;
                    double currSize = nodes[curr].size;
                    
                    double4 oldP = posMass[oldBody];
                    double3 oldPos3 = make_double3(oldP.x, oldP.y, oldP.z);
                    double oldMass = oldP.w;
                    int octOld = getOctant(currPos.x, currPos.y, currPos.z, oldPos3.x, oldPos3.y, oldPos3.z);
                    
                    double qs = currSize * 0.25;
                    double nx = currPos.x + ((octOld & 1) ? qs : -qs);
                    double ny = currPos.y + ((octOld & 2) ? qs : -qs);
                    double nz = currPos.z + ((octOld & 4) ? qs : -qs);
                    
                    int newChild = atomicAdd(nodeCounter, 1);
                    if (newChild >= maxNodes) return; 

                    initNode(nodes, newChild, make_double3(nx, ny, nz), currSize * 0.5);
                    
                    nodes[newChild].mass = oldMass;
                    nodes[newChild].center.x = oldPos3.x * oldMass;
                    nodes[newChild].center.y = oldPos3.y * oldMass;
                    nodes[newChild].center.z = oldPos3.z * oldMass;
                    nodes[newChild].bodyIdx = oldBody;
                    
                    __threadfence(); 
                    atomicCAS(&nodes[curr].children[octOld], -1, newChild);
                }
            }
            
            if (nodes[curr].bodyIdx == -2) {
                double3 currPos = nodes[curr].bMin;
                double currSize = nodes[curr].size;
                int oct = getOctant(currPos.x, currPos.y, currPos.z, pos.x, pos.y, pos.z);
                
                int child = nodes[curr].children[oct];
                if (child == -1) {
                    double qs = currSize * 0.25;
                    double nx = currPos.x + ((oct & 1) ? qs : -qs);
                    double ny = currPos.y + ((oct & 2) ? qs : -qs);
                    double nz = currPos.z + ((oct & 4) ? qs : -qs);
                    
                    int newChild = atomicAdd(nodeCounter, 1);
                    if (newChild >= maxNodes) return;

                    initNode(nodes, newChild, make_double3(nx, ny, nz), currSize * 0.5);
                    __threadfence();
                    
                    int oldChild = atomicCAS(&nodes[curr].children[oct], -1, newChild);
                    if (oldChild != -1) {
                        curr = oldChild; 
                    } else {
                        curr = newChild;
                    }
                } else {
                    curr = child;
                }
                break; 
            }
        }
    }
}

__global__ void computeForcesKernel(
    const double4* __restrict__ posMass,
    double3* __restrict__ vel,
    double3* __restrict__ acc,
    const GpuNode* __restrict__ nodes,
    int numBodies,
    double dt
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numBodies) return;

    double4 myP = posMass[idx];
    double3 myPos = make_double3(myP.x, myP.y, myP.z);
    double3 newForce = make_double3(0.0, 0.0, 0.0);

    int stack[MAX_STACK];
    int stackTop = 0;
    stack[stackTop++] = 0;

    while (stackTop > 0) {
        int nodeIdx = stack[--stackTop];
        GpuNode n = nodes[nodeIdx];

        if (n.mass <= 0.0) continue;
        double3 com = make_double3(n.center.x / n.mass, n.center.y / n.mass, n.center.z / n.mass);

        double dx = com.x - myPos.x;
        double dy = com.y - myPos.y;
        double dz = com.z - myPos.z;
        
        double distSq = dx*dx + dy*dy + dz*dz + SOFTENING_SQ; 
        double invDist = 1.0 / sqrt(distSq); 

        bool isLeaf = (n.bodyIdx >= 0); 
        
        if (isLeaf && n.bodyIdx == idx) continue;

        if (isLeaf || (n.size * invDist < THETA)) {
            double invDistCube = invDist * invDist * invDist;
            double f = (G_CONST * n.mass) * invDistCube;
            
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

void InitBarnesHutCUDA(size_t maxObjects) {
    currentCapacity = maxObjects + 4096;
    cudaMalloc(&d_posMass, currentCapacity * sizeof(double4));
    cudaMalloc(&d_vel, currentCapacity * sizeof(double3));
    cudaMalloc(&d_acc, currentCapacity * sizeof(double3));
    cudaMalloc(&d_nodes, currentCapacity * 4 * sizeof(GpuNode));
    cudaMalloc(&d_nodeCounter, sizeof(int));
    cudaMalloc(&d_minMax, 6 * sizeof(double)); 
}

void CleanupBarnesHutCUDA() {
    if (d_posMass) cudaFree(d_posMass);
    if (d_vel) cudaFree(d_vel);
    if (d_acc) cudaFree(d_acc);
    if (d_nodes) cudaFree(d_nodes);
    if (d_nodeCounter) cudaFree(d_nodeCounter);
    if (d_minMax) cudaFree(d_minMax);
}

void simulationStepBarnesHutCUDA(std::vector<Object>& objs, double dt, bool pause) {
    if (objs.empty() || pause) return;

    size_t N = objs.size();
    if (N > currentCapacity) {
        CleanupBarnesHutCUDA();
        InitBarnesHutCUDA(N * 2);
    }

    std::vector<double4> h_posMass(N);
    std::vector<double3> h_vel(N);
    
    for (size_t i = 0; i < N; ++i) {
        glm::dvec3 p = objs[i].GetPos();
        h_posMass[i] = make_double4(p.x, p.y, p.z, objs[i].mass);
        h_vel[i] = make_double3(objs[i].velocity.x, objs[i].velocity.y, objs[i].velocity.z);
    }

    cudaMemcpy(d_posMass, h_posMass.data(), N * sizeof(double4), cudaMemcpyHostToDevice);
    cudaMemcpy(d_vel, h_vel.data(), N * sizeof(double3), cudaMemcpyHostToDevice);

    int blocks = ((int)N + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

    double h_initMinMax[6] = {1e30, -1e30, 1e30, -1e30, 1e30, -1e30};
    cudaMemcpy(d_minMax, h_initMinMax, 6 * sizeof(double), cudaMemcpyHostToDevice);
    
    computeBoundingBoxKernel<<<blocks, THREADS_PER_BLOCK>>>(d_posMass, (int)N, d_minMax);
    
    double h_minMax[6];
    cudaMemcpy(h_minMax, d_minMax, 6 * sizeof(double), cudaMemcpyDeviceToHost);

    double cx = (h_minMax[0] + h_minMax[1]) * 0.5;
    double cy = (h_minMax[2] + h_minMax[3]) * 0.5;
    double cz = (h_minMax[4] + h_minMax[5]) * 0.5;
    double maxSpan = std::max({h_minMax[1] - h_minMax[0], h_minMax[3] - h_minMax[2], h_minMax[5] - h_minMax[4]}) * 1.01;

    initRootKernel<<<1, 1>>>(d_nodes, cx, cy, cz, maxSpan);

    int initialNodeCount = 1;
    cudaMemcpy(d_nodeCounter, &initialNodeCount, sizeof(int), cudaMemcpyHostToDevice);

    int maxNodes = currentCapacity * 4;
    buildTreeKernel<<<blocks, THREADS_PER_BLOCK>>>(d_posMass, d_nodes, (int)N, d_nodeCounter, maxNodes);

    computeForcesKernel<<<blocks, THREADS_PER_BLOCK>>>(d_posMass, d_vel, d_acc, d_nodes, (int)N, dt);
    integratePositionKernel<<<blocks, THREADS_PER_BLOCK>>>(d_posMass, d_vel, (int)N, dt);

    cudaMemcpy(h_posMass.data(), d_posMass, N * sizeof(double4), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_vel.data(), d_vel, N * sizeof(double3), cudaMemcpyDeviceToHost);

    for (size_t i = 0; i < N; ++i) {
        objs[i].position = glm::dvec3((double)h_posMass[i].x, (double)h_posMass[i].y, (double)h_posMass[i].z);
        objs[i].velocity = glm::dvec3((double)h_vel[i].x, (double)h_vel[i].y, (double)h_vel[i].z);
    }
}
