#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>
#include <thrust/scan.h>
#include <thrust/sort.h>

#include <algorithm>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "barnes_hut_cuda.cuh"
#include "object.hpp"
#include "physics.hpp"

#define THREADS_PER_BLOCK 256
#define MAX_STACK 64

namespace {

typedef unsigned long long morton_t;

constexpr int MORTON_BITS_PER_AXIS = 21;
constexpr int TREE_DEPTH = MORTON_BITS_PER_AXIS;
constexpr unsigned int MORTON_GRID_SIZE = 1u << 21; // 2097152
constexpr int NODE_CAPACITY_FACTOR = 16;
constexpr int MAX_FIXED_SUBSTEPS = 256;
constexpr double MAX_FIXED_STEP = 1.0 / 1024.0;
constexpr double G_CONST_D = physics::G;
constexpr int STATUS_BUILD_NODE_OVERFLOW = 1 << 0;
constexpr int STATUS_TRAVERSAL_STACK_OVERFLOW = 1 << 1;

struct __align__(16) GpuNode {
    double3 center;
    double mass;
    float3 bMin;
    float size;
    int children[8];
    int bodyBegin;
    int bodyCount;
};

static double4* d_posMass = nullptr;
static double3* d_vel = nullptr;
static double3* d_acc = nullptr;
static double3* d_newAcc = nullptr;

static GpuNode* d_nodes = nullptr;
static size_t currentNodeCapacity = 0;

static double* d_minMax = nullptr;
static morton_t* d_mortonCodes = nullptr;
static int* d_sortedBodyIndices = nullptr;
static int* d_leafFlags = nullptr;
static int* d_leafIds = nullptr;
static int* d_leafStarts = nullptr;
static int* d_leafEnds = nullptr;
static morton_t* d_leafCodes = nullptr;
static int* d_leafNodeIndices = nullptr;
static int* d_leafCount = nullptr;
static int* d_nodeCounter = nullptr;
static int* d_statusFlags = nullptr;

static size_t currentCapacity = 0;
static bool deviceStateValid = false;
static bool accInitialized = false;

__device__ void atomicMinDouble(double* addr, double value) {
    double old = *addr;
    while (old > value) {
        const double assumed = old;
        old = __longlong_as_double(
            atomicCAS(reinterpret_cast<unsigned long long*>(addr),
                      __double_as_longlong(assumed),
                      __double_as_longlong(value)));
        if (old == assumed) {
            break;
        }
    }
}

__device__ void atomicMaxDouble(double* addr, double value) {
    double old = *addr;
    while (old < value) {
        const double assumed = old;
        old = __longlong_as_double(
            atomicCAS(reinterpret_cast<unsigned long long*>(addr),
                      __double_as_longlong(assumed),
                      __double_as_longlong(value)));
        if (old == assumed) {
            break;
        }
    }
}

__device__ __forceinline__ void initNode(GpuNode* nodes, int idx, float3 cellCenter, float size) {
    nodes[idx].center = make_double3(0.0, 0.0, 0.0);
    nodes[idx].mass = 0.0;
    nodes[idx].bMin = cellCenter;
    nodes[idx].size = size;
    nodes[idx].bodyBegin = -1;
    nodes[idx].bodyCount = 0;
    #pragma unroll
    for (int i = 0; i < 8; ++i) {
        nodes[idx].children[i] = -1;
    }
}

__device__ __forceinline__ morton_t morton3D(unsigned int x, unsigned int y, unsigned int z) {
    morton_t code = 0;
    #pragma unroll
    for (int i = 0; i < 21; ++i) {
        morton_t bx = (x >> i) & 1ULL;
        morton_t by = (y >> i) & 1ULL;
        morton_t bz = (z >> i) & 1ULL;
        code |= (bx << (3 * i + 2));
        code |= (by << (3 * i + 1));
        code |= (bz << (3 * i));
    }
    return code;
}

__device__ __forceinline__ int octantFromCode(morton_t code, int level) {
    const int shift = 3 * (TREE_DEPTH - level);
    return static_cast<int>((code >> shift) & 0x7ULL);
}

__device__ __forceinline__ float3 childCellCenter(float3 center, float size, int octant) {
    const float childSize = size * 0.5f;
    const float offset = childSize * 0.5f;
    center.x += (octant & 1) ? offset : -offset;
    center.y += (octant & 2) ? offset : -offset;
    center.z += (octant & 4) ? offset : -offset;
    return center;
}

__device__ __forceinline__ float3 cellCenterFromCode(unsigned int code, float3 rootCenter, float rootSize) {
    float3 center = rootCenter;
    float size = rootSize;
    #pragma unroll
    for (int level = 1; level <= TREE_DEPTH; ++level) {
        const int oct = octantFromCode(code, level);
        center = childCellCenter(center, size, oct);
        size *= 0.5f;
    }
    return center;
}

__global__ void resetBuildStateKernel(double* dMinMax, int* leafCount, int* nodeCounter, int* statusFlags) {
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        dMinMax[0] = 1e30;
        dMinMax[1] = -1e30;
        dMinMax[2] = 1e30;
        dMinMax[3] = -1e30;
        dMinMax[4] = 1e30;
        dMinMax[5] = -1e30;
        *leafCount = 0;
        *nodeCounter = 1;
        *statusFlags = 0;
    }
}

__global__ void computeBoundingBoxKernel(const double4* posMass, int numBodies, double* dMinMax) {
    __shared__ double sMinMax[6];

    if (threadIdx.x == 0) {
        sMinMax[0] = 1e30;
        sMinMax[1] = -1e30;
        sMinMax[2] = 1e30;
        sMinMax[3] = -1e30;
        sMinMax[4] = 1e30;
        sMinMax[5] = -1e30;
    }
    __syncthreads();

    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < numBodies) {
        const double4 p = posMass[idx];
        atomicMinDouble(&sMinMax[0], p.x);
        atomicMaxDouble(&sMinMax[1], p.x);
        atomicMinDouble(&sMinMax[2], p.y);
        atomicMaxDouble(&sMinMax[3], p.y);
        atomicMinDouble(&sMinMax[4], p.z);
        atomicMaxDouble(&sMinMax[5], p.z);
    }
    __syncthreads();

    if (threadIdx.x == 0) {
        atomicMinDouble(&dMinMax[0], sMinMax[0]);
        atomicMaxDouble(&dMinMax[1], sMinMax[1]);
        atomicMinDouble(&dMinMax[2], sMinMax[2]);
        atomicMaxDouble(&dMinMax[3], sMinMax[3]);
        atomicMinDouble(&dMinMax[4], sMinMax[4]);
        atomicMaxDouble(&dMinMax[5], sMinMax[5]);
    }
}

__global__ void initRootKernel(GpuNode* nodes, const double* minMax) {
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        const double minX = minMax[0];
        const double maxX = minMax[1];
        const double minY = minMax[2];
        const double maxY = minMax[3];
        const double minZ = minMax[4];
        const double maxZ = minMax[5];

        const double spanX = maxX - minX;
        const double spanY = maxY - minY;
        const double spanZ = maxZ - minZ;
        const double maxSpan = fmax(fmax(spanX, spanY), spanZ);
        const float rootSize = static_cast<float>(fmax(maxSpan * 1.0001, 1e-6));
        const float3 center = make_float3(
            static_cast<float>((minX + maxX) * 0.5),
            static_cast<float>((minY + maxY) * 0.5),
            static_cast<float>((minZ + maxZ) * 0.5));
        initNode(nodes, 0, center, rootSize);
    }
}

__global__ void computeMortonCodesKernel(
    const double4* __restrict__ posMass,
    morton_t* __restrict__ mortonCodes,
    int* __restrict__ sortedBodyIndices,
    int numBodies,
    const double* __restrict__ minMax) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numBodies) {
        return;
    }

    const double minX = minMax[0];
    const double maxX = minMax[1];
    const double minY = minMax[2];
    const double maxY = minMax[3];
    const double minZ = minMax[4];
    const double maxZ = minMax[5];
    const double span = fmax(fmax(maxX - minX, maxY - minY), maxZ - minZ);
    const double invSpan = 1.0 / fmax(span, 1e-18);

    const double4 p = posMass[idx];
    const float nx = fminf(fmaxf(static_cast<float>((p.x - minX) * invSpan), 0.0f), 0.999999f);
    const float ny = fminf(fmaxf(static_cast<float>((p.y - minY) * invSpan), 0.0f), 0.999999f);
    const float nz = fminf(fmaxf(static_cast<float>((p.z - minZ) * invSpan), 0.0f), 0.999999f);

    const unsigned int ix = min(static_cast<unsigned int>(nx * static_cast<float>(MORTON_GRID_SIZE)), MORTON_GRID_SIZE - 1u);
    const unsigned int iy = min(static_cast<unsigned int>(ny * static_cast<float>(MORTON_GRID_SIZE)), MORTON_GRID_SIZE - 1u);
    const unsigned int iz = min(static_cast<unsigned int>(nz * static_cast<float>(MORTON_GRID_SIZE)), MORTON_GRID_SIZE - 1u);

    mortonCodes[idx] = morton3D(ix, iy, iz);
    sortedBodyIndices[idx] = idx;
}

__global__ void markLeafStartsKernel(const morton_t* mortonCodes, int* leafFlags, int numBodies) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numBodies) {
        return;
    }
    leafFlags[idx] = (idx == 0 || mortonCodes[idx] != mortonCodes[idx - 1]) ? 1 : 0;
}

__global__ void buildLeafMetadataKernel(
    const morton_t* __restrict__ mortonCodes,
    const int* __restrict__ leafFlags,
    const int* __restrict__ leafIds,
    int* __restrict__ leafStarts,
    int* __restrict__ leafEnds,
    morton_t* __restrict__ leafCodes,
    int* __restrict__ leafCount,
    int numBodies) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numBodies) {
        return;
    }

    if (leafFlags[idx]) {
        const int leafId = leafIds[idx];
        leafStarts[leafId] = idx;
        leafCodes[leafId] = mortonCodes[idx];
    }

    if (idx == numBodies - 1 || mortonCodes[idx] != mortonCodes[idx + 1]) {
        const int leafId = leafIds[idx];
        leafEnds[leafId] = idx + 1;
        if (idx == numBodies - 1) {
            *leafCount = leafId + 1;
        }
    }
}

__global__ void initLeafNodesKernel(
    const double4* __restrict__ posMass,
    const int* __restrict__ sortedBodyIndices,
    const int* __restrict__ leafStarts,
    const int* __restrict__ leafEnds,
    const morton_t* __restrict__ leafCodes,
    const int* __restrict__ leafCount,
    GpuNode* __restrict__ nodes,
    int* __restrict__ leafNodeIndices,
    int* __restrict__ nodeCounter) {
    const int leafId = blockIdx.x * blockDim.x + threadIdx.x;
    const int numLeaves = *leafCount;
    if (leafId >= numLeaves) {
        return;
    }

    const int start = leafStarts[leafId];
    const int end = leafEnds[leafId];
    const morton_t code = leafCodes[leafId];
    const int nodeIdx = leafId + 1;

    const float3 rootCenter = nodes[0].bMin;
    const float rootSize = nodes[0].size;
    const float leafSize = rootSize * (1.0f / static_cast<float>(1 << 21));
    const float3 leafCenter = cellCenterFromCode(code, rootCenter, rootSize);

    initNode(nodes, nodeIdx, leafCenter, leafSize);

    double mass = 0.0;
    double3 weighted = make_double3(0.0, 0.0, 0.0);
    for (int i = start; i < end; ++i) {
        const int bodyIdx = sortedBodyIndices[i];
        const double4 p = posMass[bodyIdx];
        const double bodyMass = p.w;
        mass += bodyMass;
        weighted.x += p.x * bodyMass;
        weighted.y += p.y * bodyMass;
        weighted.z += p.z * bodyMass;
    }

    nodes[nodeIdx].mass = mass;
    nodes[nodeIdx].center = weighted;
    nodes[nodeIdx].bodyBegin = start;
    nodes[nodeIdx].bodyCount = end - start;
    leafNodeIndices[leafId] = nodeIdx;

    if (leafId == 0) {
        *nodeCounter = numLeaves + 1;
    }
}

__global__ void buildInternalNodesKernel(
    GpuNode* __restrict__ nodes,
    const morton_t* __restrict__ leafCodes,
    const int* __restrict__ leafNodeIndices,
    const int* __restrict__ leafCount,
    int* __restrict__ nodeCounter,
    int* __restrict__ statusFlags,
    int maxNodes) {
    const int leafId = blockIdx.x * blockDim.x + threadIdx.x;
    const int numLeaves = *leafCount;
    if (leafId >= numLeaves) {
        return;
    }

    const int leafNodeIdx = leafNodeIndices[leafId];
    const morton_t code = leafCodes[leafId];
    const double leafMass = nodes[leafNodeIdx].mass;
    const double3 leafWeightedCenter = nodes[leafNodeIdx].center;

    atomicAdd(&nodes[0].mass, leafMass);
    atomicAdd(&nodes[0].center.x, leafWeightedCenter.x);
    atomicAdd(&nodes[0].center.y, leafWeightedCenter.y);
    atomicAdd(&nodes[0].center.z, leafWeightedCenter.z);

    int current = 0;
    float3 cellCenter = nodes[0].bMin;
    float cellSize = nodes[0].size;

    #pragma unroll
    for (int level = 1; level < TREE_DEPTH; ++level) {
        const int oct = octantFromCode(code, level);
        const float childSize = cellSize * 0.5f;
        const float3 nextCenter = childCellCenter(cellCenter, cellSize, oct);

        int next = nodes[current].children[oct];
        if (next == -1) {
            const int candidate = atomicAdd(nodeCounter, 1);
            if (candidate >= maxNodes) {
                atomicOr(statusFlags, STATUS_BUILD_NODE_OVERFLOW);
                return;
            }

            initNode(nodes, candidate, nextCenter, childSize);
            __threadfence();

            const int existing = atomicCAS(&nodes[current].children[oct], -1, candidate);
            next = (existing == -1) ? candidate : existing;
        }

        atomicAdd(&nodes[next].mass, leafMass);
        atomicAdd(&nodes[next].center.x, leafWeightedCenter.x);
        atomicAdd(&nodes[next].center.y, leafWeightedCenter.y);
        atomicAdd(&nodes[next].center.z, leafWeightedCenter.z);

        current = next;
        cellCenter = nextCenter;
        cellSize = childSize;
    }

    const int leafOct = octantFromCode(code, TREE_DEPTH);
    atomicCAS(&nodes[current].children[leafOct], -1, leafNodeIdx);
}

__global__ void finalizeNodesKernel(GpuNode* nodes, const int* nodeCounter) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int usedNodes = *nodeCounter;
    if (idx >= usedNodes) {
        return;
    }

    const double mass = nodes[idx].mass;
    if (mass > 0.0) {
        const double invMass = 1.0 / mass;
        nodes[idx].center.x *= invMass;
        nodes[idx].center.y *= invMass;
        nodes[idx].center.z *= invMass;
    }
}

__global__ void computeForcesKernel(
    const double4* __restrict__ posMass,
    const int* __restrict__ sortedBodyIndices,
    const GpuNode* __restrict__ nodes,
    double3* __restrict__ outAcc,
    int* __restrict__ statusFlags,
    int numBodies,
    float theta,
    double softeningSq) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numBodies) {
        return;
    }

    const double4 self = posMass[idx];
    const double3 myPos = make_double3(self.x, self.y, self.z);

    double3 acc = make_double3(0.0, 0.0, 0.0);
    int stack[MAX_STACK];
    int stackTop = 1;
    stack[0] = 0;

    while (stackTop > 0) {
        const GpuNode node = nodes[stack[--stackTop]];
        if (node.mass <= 0.0) {
            continue;
        }

        if (node.bodyCount > 0) {
            const int end = node.bodyBegin + node.bodyCount;
            for (int i = node.bodyBegin; i < end; ++i) {
                const int otherIdx = sortedBodyIndices[i];
                if (otherIdx == idx) {
                    continue;
                }

                const double4 other = posMass[otherIdx];
                const double otherMass = other.w;
                if (otherMass <= 0.0) {
                    continue;
                }

                const double dx = other.x - myPos.x;
                const double dy = other.y - myPos.y;
                const double dz = other.z - myPos.z;
                const double distSq = dx * dx + dy * dy + dz * dz + softeningSq;
                const double invDist = rsqrt(distSq);
                const double invDist3 = invDist * invDist * invDist;
                const double scale = G_CONST_D * otherMass * invDist3;

                acc.x += scale * dx;
                acc.y += scale * dy;
                acc.z += scale * dz;
            }
            continue;
        }

        const double dx = node.center.x - myPos.x;
        const double dy = node.center.y - myPos.y;
        const double dz = node.center.z - myPos.z;
        const double distSq = dx * dx + dy * dy + dz * dz + softeningSq;
        const double invDist = rsqrt(distSq);

        if (static_cast<double>(node.size) * invDist < static_cast<double>(theta)) {
            const double invDist3 = invDist * invDist * invDist;
            const double scale = G_CONST_D * node.mass * invDist3;
            acc.x += scale * dx;
            acc.y += scale * dy;
            acc.z += scale * dz;
            continue;
        }

        #pragma unroll
        for (int child = 0; child < 8; ++child) {
            const int next = node.children[child];
            if (next >= 0) {
                if (stackTop < MAX_STACK) {
                    stack[stackTop++] = next;
                } else {
                    atomicOr(statusFlags, STATUS_TRAVERSAL_STACK_OVERFLOW);
                }
            }
        }
    }

    outAcc[idx] = acc;
}

__global__ void integratePositionKernel(
    double4* posMass,
    const double3* vel,
    const double3* acc,
    int numBodies,
    double dt) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numBodies) {
        return;
    }

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
    double dt) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numBodies) {
        return;
    }

    const double halfDt = 0.5 * dt;
    vel[idx].x += (acc[idx].x + newAcc[idx].x) * halfDt;
    vel[idx].y += (acc[idx].y + newAcc[idx].y) * halfDt;
    vel[idx].z += (acc[idx].z + newAcc[idx].z) * halfDt;
    acc[idx] = newAcc[idx];
}

int fixedSubstepCount(double dt) {
    if (dt <= 0.0) {
        return 0;
    }

    const int steps = static_cast<int>(std::ceil(dt / MAX_FIXED_STEP));
    return std::max(1, std::min(MAX_FIXED_SUBSTEPS, steps));
}

void buildTreeAndComputeAcceleration(
    size_t numBodies,
    int blocks,
    float theta,
    double softeningSq,
    double3* outAcc) {
    resetBuildStateKernel<<<1, 1>>>(d_minMax, d_leafCount, d_nodeCounter, d_statusFlags);
    CHECK_CUDA(cudaGetLastError());

    computeBoundingBoxKernel<<<blocks, THREADS_PER_BLOCK>>>(d_posMass, static_cast<int>(numBodies), d_minMax);
    CHECK_CUDA(cudaGetLastError());

    initRootKernel<<<1, 1>>>(d_nodes, d_minMax);
    CHECK_CUDA(cudaGetLastError());

    computeMortonCodesKernel<<<blocks, THREADS_PER_BLOCK>>>(
        d_posMass, d_mortonCodes, d_sortedBodyIndices, static_cast<int>(numBodies), d_minMax);
    CHECK_CUDA(cudaGetLastError());

    thrust::device_ptr<morton_t> mortonBegin(d_mortonCodes);
    thrust::device_ptr<int> indexBegin(d_sortedBodyIndices);
    thrust::sort_by_key(thrust::device, mortonBegin, mortonBegin + numBodies, indexBegin);

    markLeafStartsKernel<<<blocks, THREADS_PER_BLOCK>>>(d_mortonCodes, d_leafFlags, static_cast<int>(numBodies));
    CHECK_CUDA(cudaGetLastError());

    thrust::device_ptr<int> flagBegin(d_leafFlags);
    thrust::device_ptr<int> leafIdBegin(d_leafIds);
    thrust::exclusive_scan(thrust::device, flagBegin, flagBegin + numBodies, leafIdBegin);

    CHECK_CUDA(cudaMemset(d_leafStarts, 0xFF, currentCapacity * sizeof(int)));
    CHECK_CUDA(cudaMemset(d_leafEnds, 0, currentCapacity * sizeof(int)));

    buildLeafMetadataKernel<<<blocks, THREADS_PER_BLOCK>>>(
        d_mortonCodes,
        d_leafFlags,
        d_leafIds,
        d_leafStarts,
        d_leafEnds,
        d_leafCodes,
        d_leafCount,
        static_cast<int>(numBodies));
    CHECK_CUDA(cudaGetLastError());

    initLeafNodesKernel<<<blocks, THREADS_PER_BLOCK>>>(
        d_posMass,
        d_sortedBodyIndices,
        d_leafStarts,
        d_leafEnds,
        d_leafCodes,
        d_leafCount,
        d_nodes,
        d_leafNodeIndices,
        d_nodeCounter);
    CHECK_CUDA(cudaGetLastError());

    buildInternalNodesKernel<<<blocks, THREADS_PER_BLOCK>>>(
        d_nodes,
        d_leafCodes,
        d_leafNodeIndices,
        d_leafCount,
        d_nodeCounter,
        d_statusFlags,
        static_cast<int>(currentNodeCapacity));
    CHECK_CUDA(cudaGetLastError());

    const int nodeBlocks = static_cast<int>((currentNodeCapacity + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
    finalizeNodesKernel<<<nodeBlocks, THREADS_PER_BLOCK>>>(d_nodes, d_nodeCounter);
    CHECK_CUDA(cudaGetLastError());

    computeForcesKernel<<<blocks, THREADS_PER_BLOCK>>>(
        d_posMass,
        d_sortedBodyIndices,
        d_nodes,
        outAcc,
        d_statusFlags,
        static_cast<int>(numBodies),
        theta,
        softeningSq);
    CHECK_CUDA(cudaGetLastError());

    int statusFlags = 0;
    CHECK_CUDA(cudaMemcpy(&statusFlags, d_statusFlags, sizeof(int), cudaMemcpyDeviceToHost));
    if (statusFlags != 0) {
        if ((statusFlags & STATUS_BUILD_NODE_OVERFLOW) != 0) {
            throw std::runtime_error("Barnes-Hut CUDA tree capacity exhausted during build.");
        }
        if ((statusFlags & STATUS_TRAVERSAL_STACK_OVERFLOW) != 0) {
            throw std::runtime_error("Barnes-Hut CUDA traversal stack overflowed.");
        }
    }
}

} // namespace

void InitBarnesHutCUDA(size_t maxObjects) {
    currentCapacity = maxObjects + 4096;
    currentNodeCapacity = currentCapacity * NODE_CAPACITY_FACTOR;

    CHECK_CUDA(cudaMalloc(&d_posMass, currentCapacity * sizeof(double4)));
    CHECK_CUDA(cudaMalloc(&d_vel, currentCapacity * sizeof(double3)));
    CHECK_CUDA(cudaMalloc(&d_acc, currentCapacity * sizeof(double3)));
    CHECK_CUDA(cudaMalloc(&d_newAcc, currentCapacity * sizeof(double3)));

    CHECK_CUDA(cudaMalloc(&d_nodes, currentNodeCapacity * sizeof(GpuNode)));
    CHECK_CUDA(cudaMalloc(&d_minMax, 6 * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&d_mortonCodes, currentCapacity * sizeof(morton_t)));
    CHECK_CUDA(cudaMalloc(&d_sortedBodyIndices, currentCapacity * sizeof(int)));
    CHECK_CUDA(cudaMalloc(&d_leafFlags, currentCapacity * sizeof(int)));
    CHECK_CUDA(cudaMalloc(&d_leafIds, currentCapacity * sizeof(int)));
    CHECK_CUDA(cudaMalloc(&d_leafStarts, currentCapacity * sizeof(int)));
    CHECK_CUDA(cudaMalloc(&d_leafEnds, currentCapacity * sizeof(int)));
    CHECK_CUDA(cudaMalloc(&d_leafCodes, currentCapacity * sizeof(morton_t)));
    CHECK_CUDA(cudaMalloc(&d_leafNodeIndices, currentCapacity * sizeof(int)));
    CHECK_CUDA(cudaMalloc(&d_leafCount, sizeof(int)));
    CHECK_CUDA(cudaMalloc(&d_nodeCounter, sizeof(int)));
    CHECK_CUDA(cudaMalloc(&d_statusFlags, sizeof(int)));

    deviceStateValid = false;
    accInitialized = false;
}

void CleanupBarnesHutCUDA() {
    if (d_posMass) cudaFree(d_posMass);
    if (d_vel) cudaFree(d_vel);
    if (d_acc) cudaFree(d_acc);
    if (d_newAcc) cudaFree(d_newAcc);
    if (d_nodes) cudaFree(d_nodes);
    if (d_minMax) cudaFree(d_minMax);
    if (d_mortonCodes) cudaFree(d_mortonCodes);
    if (d_sortedBodyIndices) cudaFree(d_sortedBodyIndices);
    if (d_leafFlags) cudaFree(d_leafFlags);
    if (d_leafIds) cudaFree(d_leafIds);
    if (d_leafStarts) cudaFree(d_leafStarts);
    if (d_leafEnds) cudaFree(d_leafEnds);
    if (d_leafCodes) cudaFree(d_leafCodes);
    if (d_leafNodeIndices) cudaFree(d_leafNodeIndices);
    if (d_leafCount) cudaFree(d_leafCount);
    if (d_nodeCounter) cudaFree(d_nodeCounter);
    if (d_statusFlags) cudaFree(d_statusFlags);

    d_posMass = nullptr;
    d_vel = nullptr;
    d_acc = nullptr;
    d_newAcc = nullptr;
    d_nodes = nullptr;
    d_minMax = nullptr;
    d_mortonCodes = nullptr;
    d_sortedBodyIndices = nullptr;
    d_leafFlags = nullptr;
    d_leafIds = nullptr;
    d_leafStarts = nullptr;
    d_leafEnds = nullptr;
    d_leafCodes = nullptr;
    d_leafNodeIndices = nullptr;
    d_leafCount = nullptr;
    d_nodeCounter = nullptr;
    d_statusFlags = nullptr;

    currentCapacity = 0;
    currentNodeCapacity = 0;
    deviceStateValid = false;
    accInitialized = false;
}

void simulationStepBarnesHutCUDA(std::vector<Object>& objs, double dt, bool pause, bool forceSync) {
    if (objs.empty()) {
        return;
    }
    if (pause && !forceSync) {
        return;
    }

    const size_t N = objs.size();
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
            const glm::dvec3 p = objs[i].GetPos();
            h_posMass[i] = make_double4(p.x, p.y, p.z, objs[i].mass);
            h_vel[i] = make_double3(objs[i].velocity.x, objs[i].velocity.y, objs[i].velocity.z);
            h_acc[i] = make_double3(objs[i].acceleration.x, objs[i].acceleration.y, objs[i].acceleration.z);
        }

        CHECK_CUDA(cudaMemcpy(d_posMass, h_posMass.data(), N * sizeof(double4), cudaMemcpyHostToDevice));
        CHECK_CUDA(cudaMemcpy(d_vel, h_vel.data(), N * sizeof(double3), cudaMemcpyHostToDevice));
        CHECK_CUDA(cudaMemcpy(d_acc, h_acc.data(), N * sizeof(double3), cudaMemcpyHostToDevice));
        deviceStateValid = true;
    }

    const int blocks = static_cast<int>((N + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
    const float theta = static_cast<float>(physics::getBarnesHutTheta());
    const double soft = physics::getSofteningAU();
    const double softSq = soft * soft;

    if (dt > 0.0) {
        if (!accInitialized) {
            buildTreeAndComputeAcceleration(N, blocks, theta, softSq, d_acc);
            accInitialized = true;
        }

        const int substeps = fixedSubstepCount(dt);
        const double h = dt / static_cast<double>(substeps);

        for (int step = 0; step < substeps; ++step) {
            integratePositionKernel<<<blocks, THREADS_PER_BLOCK>>>(d_posMass, d_vel, d_acc, static_cast<int>(N), h);
            CHECK_CUDA(cudaGetLastError());

            buildTreeAndComputeAcceleration(N, blocks, theta, softSq, d_newAcc);
            updateVelocityKernel<<<blocks, THREADS_PER_BLOCK>>>(d_vel, d_acc, d_newAcc, static_cast<int>(N), h);
            CHECK_CUDA(cudaGetLastError());
        }
    }

    CHECK_CUDA(cudaMemcpy(h_posMass.data(), d_posMass, N * sizeof(double4), cudaMemcpyDeviceToHost));
    for (size_t i = 0; i < N; ++i) {
        objs[i].position = glm::dvec3(h_posMass[i].x, h_posMass[i].y, h_posMass[i].z);
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
