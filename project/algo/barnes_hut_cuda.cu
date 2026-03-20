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
static constexpr double G_CONST = physics::G;

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
static double3* d_newAcc = nullptr;
static GpuNode* d_nodes = nullptr;
static int* d_nodeCounter = nullptr; 
static double* d_minMax = nullptr; 
static double* d_maxAcc = nullptr;

static size_t currentCapacity = 0;
static bool deviceStateValid = false;
static bool accInitialized = false;

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

__global__ void resetTreeStateKernel(double* d_minMax, int* d_nodeCounter) {
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        d_minMax[0] = 1e30;  d_minMax[1] = -1e30;
        d_minMax[2] = 1e30;  d_minMax[3] = -1e30;
        d_minMax[4] = 1e30;  d_minMax[5] = -1e30;
        *d_nodeCounter = 1;
    }
}

__global__ void initRootKernel(GpuNode* nodes, const double* minMax) {
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        const double minX = minMax[0], maxX = minMax[1];
        const double minY = minMax[2], maxY = minMax[3];
        const double minZ = minMax[4], maxZ = minMax[5];
        const double cx = (minX + maxX) * 0.5;
        const double cy = (minY + maxY) * 0.5;
        const double cz = (minZ + maxZ) * 0.5;
        const double maxSpan = fmax(fmax(maxX - minX, maxY - minY), maxZ - minZ) * 1.01;
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
    double3* __restrict__ outAcc,
    const GpuNode* __restrict__ nodes,
    int numBodies,
    double theta,
    double softeningSq
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
        
        double distSq = dx*dx + dy*dy + dz*dz + softeningSq; 
        if (distSq <= 0.0) continue;
        double invDist = 1.0 / sqrt(distSq); 

        bool isLeaf = (n.bodyIdx >= 0); 
        
        if (isLeaf && n.bodyIdx == idx) continue;

        if (isLeaf || (n.size * invDist < theta)) {
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

    outAcc[idx] = newForce;
}

__global__ void maxAccelKernel(const double3* acc, int n, double* outMax) {
    __shared__ double smax[THREADS_PER_BLOCK];
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    double v = 0.0;
    if (idx < n) {
        const double3 a = acc[idx];
        v = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
    }
    smax[threadIdx.x] = v;
    __syncthreads();
    for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
        if (threadIdx.x < stride) {
            smax[threadIdx.x] = fmax(smax[threadIdx.x], smax[threadIdx.x + stride]);
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        atomicMax(reinterpret_cast<unsigned long long*>(outMax), __double_as_longlong(smax[0]));
    }
}

__global__ void integratePositionKernel(
    double4* posMass,
    const double3* vel,
    const double3* acc,
    int numBodies,
    double dt
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numBodies) return;

    const double halfDt2 = 0.5 * dt * dt;
    posMass[idx].x += vel[idx].x * dt + acc[idx].x * halfDt2;
    posMass[idx].y += vel[idx].y * dt + acc[idx].y * halfDt2;
    posMass[idx].z += vel[idx].z * dt + acc[idx].z * halfDt2;
}

__global__ void updateVelocityKernel(
    double3* vel,
    double3* acc,
    const double3* newAcc,
    int numBodies,
    double dt
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numBodies) return;

    const double halfDt = 0.5 * dt;
    vel[idx].x += (acc[idx].x + newAcc[idx].x) * halfDt;
    vel[idx].y += (acc[idx].y + newAcc[idx].y) * halfDt;
    vel[idx].z += (acc[idx].z + newAcc[idx].z) * halfDt;
    acc[idx] = newAcc[idx];
}

void InitBarnesHutCUDA(size_t maxObjects) {
    currentCapacity = maxObjects + 4096;
    CHECK_CUDA(cudaMalloc(&d_posMass, currentCapacity * sizeof(double4)));
    CHECK_CUDA(cudaMalloc(&d_vel, currentCapacity * sizeof(double3)));
    CHECK_CUDA(cudaMalloc(&d_acc, currentCapacity * sizeof(double3)));
    CHECK_CUDA(cudaMalloc(&d_newAcc, currentCapacity * sizeof(double3)));
    CHECK_CUDA(cudaMalloc(&d_nodes, currentCapacity * 4 * sizeof(GpuNode)));
    CHECK_CUDA(cudaMalloc(&d_nodeCounter, sizeof(int)));
    CHECK_CUDA(cudaMalloc(&d_minMax, 6 * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&d_maxAcc, sizeof(double)));
    deviceStateValid = false;
    accInitialized = false;
}

void CleanupBarnesHutCUDA() {
    if (d_posMass) cudaFree(d_posMass);
    if (d_vel) cudaFree(d_vel);
    if (d_acc) cudaFree(d_acc);
    if (d_newAcc) cudaFree(d_newAcc);
    if (d_nodes) cudaFree(d_nodes);
    if (d_nodeCounter) cudaFree(d_nodeCounter);
    if (d_minMax) cudaFree(d_minMax);
    if (d_maxAcc) cudaFree(d_maxAcc);
    d_posMass = nullptr;
    d_vel = nullptr;
    d_acc = nullptr;
    d_newAcc = nullptr;
    d_nodes = nullptr;
    d_nodeCounter = nullptr;
    d_minMax = nullptr;
    d_maxAcc = nullptr;
    currentCapacity = 0;
    deviceStateValid = false;
    accInitialized = false;
}

void simulationStepBarnesHutCUDA(std::vector<Object>& objs, double dt, bool pause, bool forceSync) {
    if (objs.empty()) return;
    if (pause && !forceSync) return;

    size_t N = objs.size();
    if (N > currentCapacity) {
        CleanupBarnesHutCUDA();
        InitBarnesHutCUDA(N * 2);
    }

    static std::vector<double4> h_posMass;
    static std::vector<double3> h_vel;
    static std::vector<double3> h_acc;
    h_posMass.resize(N);
    h_vel.resize(N);
    h_acc.resize(N);
    
    if (!deviceStateValid) {
        for (size_t i = 0; i < N; ++i) {
            glm::dvec3 p = objs[i].GetPos();
            h_posMass[i] = make_double4(p.x, p.y, p.z, objs[i].mass);
            h_vel[i] = make_double3(objs[i].velocity.x, objs[i].velocity.y, objs[i].velocity.z);
            h_acc[i] = make_double3(objs[i].acceleration.x, objs[i].acceleration.y, objs[i].acceleration.z);
        }

        CHECK_CUDA(cudaMemcpy(d_posMass, h_posMass.data(), N * sizeof(double4), cudaMemcpyHostToDevice));
        CHECK_CUDA(cudaMemcpy(d_vel, h_vel.data(), N * sizeof(double3), cudaMemcpyHostToDevice));
        CHECK_CUDA(cudaMemcpy(d_acc, h_acc.data(), N * sizeof(double3), cudaMemcpyHostToDevice));
        deviceStateValid = true;
    }

    int blocks = ((int)N + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    const double theta = physics::getBarnesHutTheta();
    const double soft = physics::getSofteningAU();
    const double softSq = soft * soft;

    auto buildTree = [&](double4* positions) {
        resetTreeStateKernel<<<1, 1>>>(d_minMax, d_nodeCounter);
        CHECK_CUDA(cudaGetLastError());
        computeBoundingBoxKernel<<<blocks, THREADS_PER_BLOCK>>>(positions, (int)N, d_minMax);
        CHECK_CUDA(cudaGetLastError());
        initRootKernel<<<1, 1>>>(d_nodes, d_minMax);
        CHECK_CUDA(cudaGetLastError());

        int maxNodes = static_cast<int>(currentCapacity * 4);
        buildTreeKernel<<<blocks, THREADS_PER_BLOCK>>>(positions, d_nodes, (int)N, d_nodeCounter, maxNodes);
        CHECK_CUDA(cudaGetLastError());
    };

    if (dt > 0.0) {
        double remaining = dt;
        int substeps = 0;
        constexpr int MAX_SUBSTEPS = 1024;
        while (remaining > 0.0 && substeps < MAX_SUBSTEPS) {
            // Инициализация ускорений на первом шаге для корректного Verlet.
            if (!accInitialized) {
                buildTree(d_posMass);
                computeForcesKernel<<<blocks, THREADS_PER_BLOCK>>>(d_posMass, d_acc, d_nodes, (int)N, theta, softSq);
                CHECK_CUDA(cudaGetLastError());
                accInitialized = true;
            }

            double zero = 0.0;
            CHECK_CUDA(cudaMemcpy(d_maxAcc, &zero, sizeof(double), cudaMemcpyHostToDevice));
            maxAccelKernel<<<blocks, THREADS_PER_BLOCK>>>(d_acc, (int)N, d_maxAcc);
            CHECK_CUDA(cudaGetLastError());
            double maxAcc = 0.0;
            CHECK_CUDA(cudaMemcpy(&maxAcc, d_maxAcc, sizeof(double), cudaMemcpyDeviceToHost));

            double suggested = remaining;
            if (maxAcc > 0.0) {
                suggested = std::min(suggested, 0.02 / std::sqrt(maxAcc));
            }
            double minStep = dt / 4096.0;
            double h = std::clamp(suggested, minStep, remaining);

            integratePositionKernel<<<blocks, THREADS_PER_BLOCK>>>(d_posMass, d_vel, d_acc, (int)N, h);
            CHECK_CUDA(cudaGetLastError());

            buildTree(d_posMass);
            computeForcesKernel<<<blocks, THREADS_PER_BLOCK>>>(d_posMass, d_newAcc, d_nodes, (int)N, theta, softSq);
            CHECK_CUDA(cudaGetLastError());

            updateVelocityKernel<<<blocks, THREADS_PER_BLOCK>>>(d_vel, d_acc, d_newAcc, (int)N, h);
            CHECK_CUDA(cudaGetLastError());

            remaining -= h;
            ++substeps;
        }
    }

    CHECK_CUDA(cudaMemcpy(h_posMass.data(), d_posMass, N * sizeof(double4), cudaMemcpyDeviceToHost));

    for (size_t i = 0; i < N; ++i) {
        objs[i].position = glm::dvec3((double)h_posMass[i].x, (double)h_posMass[i].y, (double)h_posMass[i].z);
    }

    if (forceSync) {
        CHECK_CUDA(cudaMemcpy(h_vel.data(), d_vel, N * sizeof(double3), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(h_acc.data(), d_acc, N * sizeof(double3), cudaMemcpyDeviceToHost));
        for (size_t i = 0; i < N; ++i) {
            objs[i].velocity = glm::dvec3(h_vel[i].x, h_vel[i].y, h_vel[i].z);
            objs[i].acceleration = glm::dvec3(h_acc[i].x, h_acc[i].y, h_acc[i].z);
        }
    }
}
