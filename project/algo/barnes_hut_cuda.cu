#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include "barnes_hut_cuda.cuh"
#include "object.hpp"

#define THREADS_PER_BLOCK 256
#define MAX_STACK 64
#define THETA 0.5f
#define G_CONST 6.6743e-20f
#define SOFTENING_SQ 0.1f 

struct __align__(16) GpuNode {
    float3 center; 
    float mass;
    float3 bMin;   
    float size;
    int children[8];
    int bodyIdx;
};

static float4* d_posMass = nullptr;
static float3* d_vel = nullptr;
static float3* d_acc = nullptr;
static GpuNode* d_nodes = nullptr;
static int* d_nodeCounter = nullptr; 
static float* d_minMax = nullptr; 

static size_t currentCapacity = 0;

__device__ void atomicMinFloat(float* addr, float value) {
    float old = *addr, assumed;
    if (old <= value) return;
    do {
        assumed = old;
        old = __int_as_float(atomicCAS((int*)addr, __float_as_int(assumed), __float_as_int(value)));
    } while (assumed != old && old > value);
}

__device__ void atomicMaxFloat(float* addr, float value) {
    float old = *addr, assumed;
    if (old >= value) return;
    do {
        assumed = old;
        old = __int_as_float(atomicCAS((int*)addr, __float_as_int(assumed), __float_as_int(value)));
    } while (assumed != old && old < value);
}

__device__ int getOctant(float tx, float ty, float tz, float px, float py, float pz) {
    int idx = 0;
    if (px >= tx) idx |= 1;
    if (py >= ty) idx |= 2;
    if (pz >= tz) idx |= 4;
    return idx;
}

__device__ void initNode(GpuNode* nodes, int idx, float3 bMin, float size) {
    nodes[idx].bMin = bMin;
    nodes[idx].size = size;
    nodes[idx].mass = 0.0f;
    nodes[idx].center = make_float3(0.0f, 0.0f, 0.0f);
    nodes[idx].bodyIdx = -1;
    #pragma unroll
    for (int i = 0; i < 8; ++i) {
        nodes[idx].children[i] = -1;
    }
}

__global__ void computeBoundingBoxKernel(const float4* posMass, int numBodies, float* d_minMax) {
    __shared__ float s_minMax[6]; 

    if (threadIdx.x == 0) {
        s_minMax[0] = 1e30f; s_minMax[1] = -1e30f;
        s_minMax[2] = 1e30f; s_minMax[3] = -1e30f;
        s_minMax[4] = 1e30f; s_minMax[5] = -1e30f;
    }
    __syncthreads();

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < numBodies) {
        float4 p = posMass[idx];
        atomicMinFloat(&s_minMax[0], p.x); atomicMaxFloat(&s_minMax[1], p.x);
        atomicMinFloat(&s_minMax[2], p.y); atomicMaxFloat(&s_minMax[3], p.y);
        atomicMinFloat(&s_minMax[4], p.z); atomicMaxFloat(&s_minMax[5], p.z);
    }
    __syncthreads();

    if (threadIdx.x == 0) {
        atomicMinFloat(&d_minMax[0], s_minMax[0]); atomicMaxFloat(&d_minMax[1], s_minMax[1]);
        atomicMinFloat(&d_minMax[2], s_minMax[2]); atomicMaxFloat(&d_minMax[3], s_minMax[3]);
        atomicMinFloat(&d_minMax[4], s_minMax[4]); atomicMaxFloat(&d_minMax[5], s_minMax[5]);
    }
}

__global__ void initRootKernel(GpuNode* nodes, float cx, float cy, float cz, float maxSpan) {
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        initNode(nodes, 0, make_float3(cx, cy, cz), maxSpan);
    }
}

__global__ void buildTreeKernel(const float4* __restrict__ posMass, GpuNode* __restrict__ nodes, int numBodies, int* nodeCounter, int maxNodes) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numBodies) return;

    float4 p = posMass[i];
    float3 pos = make_float3(p.x, p.y, p.z);
    float mass = p.w;
    if (mass <= 0.0f) return;

    int curr = 0; 
    
    while (true) {
        atomicAdd(&nodes[curr].mass, mass);
        atomicAdd(&nodes[curr].center.x, pos.x * mass);
        atomicAdd(&nodes[curr].center.y, pos.y * mass);
        atomicAdd(&nodes[curr].center.z, pos.z * mass);
        
        if (nodes[curr].size < 1e-5f) return; 
        
        while (true) {
            int oldBody = atomicCAS(&nodes[curr].bodyIdx, -1, i);
            
            if (oldBody == -1 || oldBody == i) return; 
            
            if (oldBody >= 0) {
                if (atomicCAS(&nodes[curr].bodyIdx, oldBody, -2) == oldBody) {
                    float3 currPos = nodes[curr].bMin;
                    float currSize = nodes[curr].size;
                    
                    float4 oldP = posMass[oldBody];
                    float3 oldPos3 = make_float3(oldP.x, oldP.y, oldP.z);
                    float oldMass = oldP.w;
                    int octOld = getOctant(currPos.x, currPos.y, currPos.z, oldPos3.x, oldPos3.y, oldPos3.z);
                    
                    float qs = currSize * 0.25f;
                    float nx = currPos.x + ((octOld & 1) ? qs : -qs);
                    float ny = currPos.y + ((octOld & 2) ? qs : -qs);
                    float nz = currPos.z + ((octOld & 4) ? qs : -qs);
                    
                    int newChild = atomicAdd(nodeCounter, 1);
                    if (newChild >= maxNodes) return; 

                    initNode(nodes, newChild, make_float3(nx, ny, nz), currSize * 0.5f);
                    
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
                float3 currPos = nodes[curr].bMin;
                float currSize = nodes[curr].size;
                int oct = getOctant(currPos.x, currPos.y, currPos.z, pos.x, pos.y, pos.z);
                
                int child = nodes[curr].children[oct];
                if (child == -1) {
                    float qs = currSize * 0.25f;
                    float nx = currPos.x + ((oct & 1) ? qs : -qs);
                    float ny = currPos.y + ((oct & 2) ? qs : -qs);
                    float nz = currPos.z + ((oct & 4) ? qs : -qs);
                    
                    int newChild = atomicAdd(nodeCounter, 1);
                    if (newChild >= maxNodes) return;

                    initNode(nodes, newChild, make_float3(nx, ny, nz), currSize * 0.5f);
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
    const float4* __restrict__ posMass,
    float3* __restrict__ vel,
    float3* __restrict__ acc,
    const GpuNode* __restrict__ nodes,
    int numBodies,
    float dt
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numBodies) return;

    float4 myP = posMass[idx];
    float3 myPos = make_float3(myP.x, myP.y, myP.z);
    float3 newForce = make_float3(0.0f, 0.0f, 0.0f);

    int stack[MAX_STACK];
    int stackTop = 0;
    stack[stackTop++] = 0;

    while (stackTop > 0) {
        int nodeIdx = stack[--stackTop];
        GpuNode n = nodes[nodeIdx];

        float3 com = make_float3(n.center.x / n.mass, n.center.y / n.mass, n.center.z / n.mass);

        float dx = com.x - myPos.x;
        float dy = com.y - myPos.y;
        float dz = com.z - myPos.z;
        
        float distSq = dx*dx + dy*dy + dz*dz + SOFTENING_SQ; 
        float invDist = rsqrtf(distSq); 

        bool isLeaf = (n.bodyIdx >= 0); 
        
        if (isLeaf && n.bodyIdx == idx) continue;

        if (isLeaf || (n.size * invDist < THETA)) {
            float invDistCube = invDist * invDist * invDist;
            float f = (G_CONST * n.mass) * invDistCube;
            
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
    float4* posMass,
    const float3* vel,
    int numBodies,
    float dt
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numBodies) return;

    posMass[idx].x += vel[idx].x * dt;
    posMass[idx].y += vel[idx].y * dt;
    posMass[idx].z += vel[idx].z * dt;
}

void InitBarnesHutCUDA(size_t maxObjects) {
    currentCapacity = maxObjects + 4096;
    cudaMalloc(&d_posMass, currentCapacity * sizeof(float4));
    cudaMalloc(&d_vel, currentCapacity * sizeof(float3));
    cudaMalloc(&d_acc, currentCapacity * sizeof(float3));
    cudaMalloc(&d_nodes, currentCapacity * 4 * sizeof(GpuNode));
    cudaMalloc(&d_nodeCounter, sizeof(int));
    cudaMalloc(&d_minMax, 6 * sizeof(float)); 
}

void CleanupBarnesHutCUDA() {
    if (d_posMass) cudaFree(d_posMass);
    if (d_vel) cudaFree(d_vel);
    if (d_acc) cudaFree(d_acc);
    if (d_nodes) cudaFree(d_nodes);
    if (d_nodeCounter) cudaFree(d_nodeCounter);
    if (d_minMax) cudaFree(d_minMax);
}

void simulationStepBarnesHutCUDA(std::vector<Object>& objs, float dt, bool pause) {
    if (objs.empty() || pause) return;

    size_t N = objs.size();
    if (N > currentCapacity) {
        CleanupBarnesHutCUDA();
        InitBarnesHutCUDA(N * 2);
    }

    std::vector<float4> h_posMass(N);
    std::vector<float3> h_vel(N);
    
    for (size_t i = 0; i < N; ++i) {
        glm::dvec3 p = objs[i].GetPos();
        h_posMass[i] = make_float4((float)p.x, (float)p.y, (float)p.z, (float)objs[i].mass);
        h_vel[i] = make_float3((float)objs[i].velocity.x, (float)objs[i].velocity.y, (float)objs[i].velocity.z);
    }

    cudaMemcpy(d_posMass, h_posMass.data(), N * sizeof(float4), cudaMemcpyHostToDevice);
    cudaMemcpy(d_vel, h_vel.data(), N * sizeof(float3), cudaMemcpyHostToDevice);

    int blocks = ((int)N + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

    float h_initMinMax[6] = {1e30f, -1e30f, 1e30f, -1e30f, 1e30f, -1e30f};
    cudaMemcpy(d_minMax, h_initMinMax, 6 * sizeof(float), cudaMemcpyHostToDevice);
    
    computeBoundingBoxKernel<<<blocks, THREADS_PER_BLOCK>>>(d_posMass, (int)N, d_minMax);
    
    float h_minMax[6];
    cudaMemcpy(h_minMax, d_minMax, 6 * sizeof(float), cudaMemcpyDeviceToHost);

    float cx = (h_minMax[0] + h_minMax[1]) * 0.5f;
    float cy = (h_minMax[2] + h_minMax[3]) * 0.5f;
    float cz = (h_minMax[4] + h_minMax[5]) * 0.5f;
    float maxSpan = std::max({h_minMax[1] - h_minMax[0], h_minMax[3] - h_minMax[2], h_minMax[5] - h_minMax[4]}) * 1.01f;

    initRootKernel<<<1, 1>>>(d_nodes, cx, cy, cz, maxSpan);

    int initialNodeCount = 1;
    cudaMemcpy(d_nodeCounter, &initialNodeCount, sizeof(int), cudaMemcpyHostToDevice);

    int maxNodes = currentCapacity * 4;
    buildTreeKernel<<<blocks, THREADS_PER_BLOCK>>>(d_posMass, d_nodes, (int)N, d_nodeCounter, maxNodes);

    computeForcesKernel<<<blocks, THREADS_PER_BLOCK>>>(d_posMass, d_vel, d_acc, d_nodes, (int)N, dt);
    integratePositionKernel<<<blocks, THREADS_PER_BLOCK>>>(d_posMass, d_vel, (int)N, dt);

    cudaMemcpy(h_posMass.data(), d_posMass, N * sizeof(float4), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_vel.data(), d_vel, N * sizeof(float3), cudaMemcpyDeviceToHost);

    for (size_t i = 0; i < N; ++i) {
        objs[i].position = glm::dvec3((double)h_posMass[i].x, (double)h_posMass[i].y, (double)h_posMass[i].z);
        objs[i].velocity = glm::dvec3((double)h_vel[i].x, (double)h_vel[i].y, (double)h_vel[i].z);
    }
}