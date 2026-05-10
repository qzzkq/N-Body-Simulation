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
#define WARPS_PER_BLOCK (THREADS_PER_BLOCK / 32)
#define WARP_STACK_DEPTH 128     // 21 уровень × 7 siblings ≈ 147 worst-case;
                                 // в пракатике с Morton-локальностью warp'а
                                 // глубина < 64. Overflow ловится через STATUS_TRAVERSAL_STACK_OVERFLOW.
#define ALWAYS_DIRECT_FACTOR_SQ 9.0f  // открывать узел если r² < этого × softening²

namespace {

typedef unsigned long long morton_t;

constexpr int MORTON_BITS_PER_AXIS = 21;
constexpr int TREE_DEPTH = MORTON_BITS_PER_AXIS;
constexpr unsigned int MORTON_GRID_SIZE = 1u << 21; // 2097152
constexpr int NODE_CAPACITY_FACTOR = 16;
// Adaptive timestep:
//   h = clamp(eta/sqrt(max|a|), dt/cap, remaining)
// Параметры eta и cap берутся в runtime из physics::getAdaptiveEta() /
// physics::getSubstepCap() — задаются через --adaptive-eta / --substep-cap
// в CLI или интерактивный prompt «Качество». Дефолт = balanced пресет.
constexpr float G_CONST_F = static_cast<float>(physics::G);
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
    // Traceless symmetric quadrupole (Qzz = -(Qxx+Qyy), хранится 5 компонент).
    // Накапливается AtomicAdd'ом из листьев со Steiner-сдвигом до COM узла.
    float Qxx;
    float Qyy;
    float Qxy;
    float Qxz;
    float Qyz;
    // Dehnen-δ: максимум ||body - COM_node|| по всем телам в поддереве.
    // MAC: открыть узел если (r - δ) ≤ size/θ. На листьях считается напрямую,
    // на внутренних узлах накапливается atomicMax-ом сверху от листьев.
    float delta;
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
static double* d_maxAcc = nullptr;
// Pinned host-side приёмник для maxAcc — убирает накладные расходы на pageable
// memory при D→H transfer внутри substep-петли. Аллоцируется через cudaMallocHost.
static double* h_maxAccPinned = nullptr;
// Pinned host-side приёмник для status flags — аналогично, читается раз в outer-step.
static int* h_statusFlagsPinned = nullptr;

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

// Block-reduce |a| → atomicMax into outMax. Bit-trick atomicMax on unsigned long long
// works for non-negative doubles since IEEE-754 ordering matches integer ordering.
__global__ void maxAccelKernel(const double3* __restrict__ acc, int n, double* outMax) {
    __shared__ double smax[THREADS_PER_BLOCK];
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
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
        atomicMax(reinterpret_cast<unsigned long long*>(outMax),
                  static_cast<unsigned long long>(__double_as_longlong(smax[0])));
    }
}

__device__ __forceinline__ void initNode(GpuNode* nodes, int idx, float3 cellCenter, float size) {
    nodes[idx].center = make_double3(0.0, 0.0, 0.0);
    nodes[idx].mass = 0.0;
    nodes[idx].bMin = cellCenter;
    nodes[idx].size = size;
    nodes[idx].bodyBegin = -1;
    nodes[idx].bodyCount = 0;
    nodes[idx].Qxx = 0.0f;
    nodes[idx].Qyy = 0.0f;
    nodes[idx].Qxy = 0.0f;
    nodes[idx].Qxz = 0.0f;
    nodes[idx].Qyz = 0.0f;
    nodes[idx].delta = 0.0f;
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

__global__ void resetBuildStateKernel(double* dMinMax, int* leafCount, int* nodeCounter) {
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        dMinMax[0] = 1e30;
        dMinMax[1] = -1e30;
        dMinMax[2] = 1e30;
        dMinMax[3] = -1e30;
        dMinMax[4] = 1e30;
        dMinMax[5] = -1e30;
        *leafCount = 0;
        *nodeCounter = 1;
        // statusFlags НЕ обнуляем здесь — это делает host один раз на outer-step,
        // а в buildTreeAndComputeAcceleration флаги атомарно ORятся между substep'ами.
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

// Считает квадруполь листа (узла с bodyCount > 0) относительно COM этого листа.
// Должен запускаться ПОСЛЕ finalizeNodesKernel — нужен корректный COM.
__global__ void computeLeafQuadrupoleKernel(
    GpuNode* __restrict__ nodes,
    const double4* __restrict__ posMass,
    const int* __restrict__ sortedBodyIndices,
    const int* __restrict__ leafNodeIndices,
    const int* __restrict__ leafCount) {
    const int leafId = blockIdx.x * blockDim.x + threadIdx.x;
    const int numLeaves = *leafCount;
    if (leafId >= numLeaves) return;

    const int nodeIdx = leafNodeIndices[leafId];
    const double3 com  = nodes[nodeIdx].center;
    const int begin = nodes[nodeIdx].bodyBegin;
    const int end   = begin + nodes[nodeIdx].bodyCount;

    double Qxx = 0.0, Qyy = 0.0, Qxy = 0.0, Qxz = 0.0, Qyz = 0.0;
    double maxR2 = 0.0;
    for (int i = begin; i < end; ++i) {
        const int bodyIdx = sortedBodyIndices[i];
        const double4 p = posMass[bodyIdx];
        const double m = p.w;
        const double dx = p.x - com.x;
        const double dy = p.y - com.y;
        const double dz = p.z - com.z;
        const double d2 = dx*dx + dy*dy + dz*dz;
        if (d2 > maxR2) maxR2 = d2;
        Qxx += m * (3.0*dx*dx - d2);
        Qyy += m * (3.0*dy*dy - d2);
        Qxy += m * (3.0*dx*dy);
        Qxz += m * (3.0*dx*dz);
        Qyz += m * (3.0*dy*dz);
    }
    nodes[nodeIdx].Qxx = static_cast<float>(Qxx);
    nodes[nodeIdx].Qyy = static_cast<float>(Qyy);
    nodes[nodeIdx].Qxy = static_cast<float>(Qxy);
    nodes[nodeIdx].Qxz = static_cast<float>(Qxz);
    nodes[nodeIdx].Qyz = static_cast<float>(Qyz);
    nodes[nodeIdx].delta = static_cast<float>(sqrt(maxR2));
}

// Накопление квадруполя вверх по дереву: каждый лист идёт от корня к себе
// и atomicAdd'ом добавляет в каждый ancestor свой вклад со Steiner-сдвигом
// (parallel-axis theorem для traceless квадруполя):
//   Q_anc[ij] += Q_leaf[ij] + m_leaf · (3 d^i d^j - δ^ij |d|²),  d = COM_leaf - COM_anc.
__global__ void accumulateQuadrupoleKernel(
    GpuNode* __restrict__ nodes,
    const morton_t* __restrict__ leafCodes,
    const int* __restrict__ leafNodeIndices,
    const int* __restrict__ leafCount) {
    const int leafId = blockIdx.x * blockDim.x + threadIdx.x;
    const int numLeaves = *leafCount;
    if (leafId >= numLeaves) return;

    const int leafNodeIdx = leafNodeIndices[leafId];
    const double3 comL    = nodes[leafNodeIdx].center;
    const double  mL      = nodes[leafNodeIdx].mass;
    const float QLxx = nodes[leafNodeIdx].Qxx;
    const float QLyy = nodes[leafNodeIdx].Qyy;
    const float QLxy = nodes[leafNodeIdx].Qxy;
    const float QLxz = nodes[leafNodeIdx].Qxz;
    const float QLyz = nodes[leafNodeIdx].Qyz;
    const float deltaL = nodes[leafNodeIdx].delta;

    const morton_t code = leafCodes[leafId];

    int current = 0;
    for (int level = 1; level <= TREE_DEPTH; ++level) {
        // Steiner shift до COM текущего ancestor'а.
        const double3 comA = nodes[current].center;
        const double dx = comL.x - comA.x;
        const double dy = comL.y - comA.y;
        const double dz = comL.z - comA.z;
        const double d2 = dx*dx + dy*dy + dz*dz;
        const float addXx = QLxx + static_cast<float>(mL * (3.0*dx*dx - d2));
        const float addYy = QLyy + static_cast<float>(mL * (3.0*dy*dy - d2));
        const float addXy = QLxy + static_cast<float>(mL * (3.0*dx*dy));
        const float addXz = QLxz + static_cast<float>(mL * (3.0*dx*dz));
        const float addYz = QLyz + static_cast<float>(mL * (3.0*dy*dz));
        atomicAdd(&nodes[current].Qxx, addXx);
        atomicAdd(&nodes[current].Qyy, addYy);
        atomicAdd(&nodes[current].Qxy, addXy);
        atomicAdd(&nodes[current].Qxz, addXz);
        atomicAdd(&nodes[current].Qyz, addYz);

        // Dehnen-δ: max ||body - COM_anc|| ≥ ||COM_leaf - COM_anc|| + δ_leaf.
        // atomicMax на float через unsigned int — корректно для неотрицательных значений.
        const float candidateDelta = static_cast<float>(sqrt(d2)) + deltaL;
        atomicMax(reinterpret_cast<unsigned int*>(&nodes[current].delta),
                  __float_as_uint(candidateDelta));

        const int oct = octantFromCode(code, level);
        const int next = nodes[current].children[oct];
        if (next < 0 || next == leafNodeIdx) break;
        current = next;
    }
}

// Warp-collaborative traversal (Burtscher-Pingali 2011) с mixed precision.
//
// Архитектура:
//   - 1 поток = 1 тело, как и раньше.
//   - Стек вынесен в shared memory: один на варп (32 потока), не на поток.
//   - Узел грузится единожды на варп — 32 потока читают одну ту же ячейку,
//     hardware-coalescing превращает это в одну транзакцию.
//   - Решение «принять monopole vs открыть»: __ballot_sync. Если хоть один
//     поток в варпе хочет раскрыть узел — варп идёт в детей и НИКТО не
//     накапливает на этом уровне (стандарт BP — стрикт-MAC по варпу).
//
// Точность:
//   - FP64 вычитание дельт (cancellation-safe для координат до 10⁵ AU),
//     далее FP32 + hardware rsqrtf, аккумулятор FP64.
//   - Always-direct: узел никогда не аппроксимируется монополем при
//     r² < ALWAYS_DIRECT_FACTOR_SQ · softening² — защита от резких ошибок
//     в плотных сжатиях (типа коллизии ядер).
__global__ void computeForcesKernel(
    const double4* __restrict__ posMass,
    const int* __restrict__ sortedBodyIndices,
    const GpuNode* __restrict__ nodes,
    double3* __restrict__ outAcc,
    int* __restrict__ statusFlags,
    int numBodies,
    float theta,
    double softeningSq) {
    __shared__ int sStack[WARPS_PER_BLOCK * WARP_STACK_DEPTH];
    __shared__ int sStackTop[WARPS_PER_BLOCK];

    const int tid     = threadIdx.x;
    const int warpId  = tid >> 5;
    const int laneId  = tid & 31;
    const int sortPos = blockIdx.x * blockDim.x + tid;

    int* myStack = &sStack[warpId * WARP_STACK_DEPTH];

    // Идём в Morton-порядке: соседние lane'ы в варпе = соседние в пространстве.
    // Это даёт высокий warp-consensus по MAC и резко снижает число открытий узла.
    const bool active = (sortPos < numBodies);
    int idx = -1;
    double myPosX = 0.0, myPosY = 0.0, myPosZ = 0.0;
    if (active) {
        idx = sortedBodyIndices[sortPos];
        const double4 self = posMass[idx];
        myPosX = self.x;
        myPosY = self.y;
        myPosZ = self.z;
    }

    const float thetaSq = theta * theta;
    const float softSqF = static_cast<float>(softeningSq);
    const float alwaysDirectR2 = softSqF * ALWAYS_DIRECT_FACTOR_SQ;

    double3 acc = make_double3(0.0, 0.0, 0.0);

    if (laneId == 0) {
        myStack[0] = 0;
        sStackTop[warpId] = 1;
    }
    __syncwarp();

    while (sStackTop[warpId] > 0) {
        int nodeIdx = -1;
        if (laneId == 0) {
            nodeIdx = myStack[--sStackTop[warpId]];
        }
        nodeIdx = __shfl_sync(0xffffffff, nodeIdx, 0);

        // Все 32 потока читают одну и ту же ячейку — hardware кэширует это
        // как одну широковещательную транзакцию.
        const GpuNode node = nodes[nodeIdx];
        if (node.mass <= 0.0) {
            continue;
        }

        // Лист: warp-cooperative — лейны блоками по 32 совместно тянут
        // данные «других» тел, потом каждый поток применяет силу к СВОЕМУ телу.
        // ВАЖНО: все 32 лейна должны вызывать __shfl_sync, поэтому проверка
        // `active` применяется только к НАКОПЛЕНИЮ, а не к самим shfl-операциям.
        if (node.bodyCount > 0) {
            const int begin = node.bodyBegin;
            const int end   = begin + node.bodyCount;
            for (int base = begin; base < end; base += 32) {
                const int slot = base + laneId;
                int otherIdx_l = -1;
                double oXl = 0.0, oYl = 0.0, oZl = 0.0, oMl = 0.0;
                if (slot < end) {
                    otherIdx_l = sortedBodyIndices[slot];
                    const double4 p = posMass[otherIdx_l];
                    oXl = p.x; oYl = p.y; oZl = p.z; oMl = p.w;
                }
                const int tile = min(32, end - base);
                #pragma unroll 8
                for (int t = 0; t < tile; ++t) {
                    // Все 32 лейна участвуют в shfl независимо от `active`.
                    const int    oIdx  = __shfl_sync(0xffffffff, otherIdx_l, t);
                    const double oX    = __shfl_sync(0xffffffff, oXl, t);
                    const double oY    = __shfl_sync(0xffffffff, oYl, t);
                    const double oZ    = __shfl_sync(0xffffffff, oZl, t);
                    const double oMass = __shfl_sync(0xffffffff, oMl, t);
                    if (!active)            continue;
                    if (oIdx == idx)        continue;
                    if (oMass <= 0.0)       continue;
                    const float dx = static_cast<float>(oX - myPosX);
                    const float dy = static_cast<float>(oY - myPosY);
                    const float dz = static_cast<float>(oZ - myPosZ);
                    const float r2 = dx * dx + dy * dy + dz * dz + softSqF;
                    const float invR  = rsqrtf(r2);
                    const float invR3 = invR * invR * invR;
                    const float scale = G_CONST_F * static_cast<float>(oMass) * invR3;
                    acc.x += static_cast<double>(scale * dx);
                    acc.y += static_cast<double>(scale * dy);
                    acc.z += static_cast<double>(scale * dz);
                }
            }
            continue;
        }

        float r2_thread = 1.0f;
        float dxF = 0.0f, dyF = 0.0f, dzF = 0.0f;
        bool acceptMe = true;
        if (active) {
            dxF = static_cast<float>(node.center.x - myPosX);
            dyF = static_cast<float>(node.center.y - myPosY);
            dzF = static_cast<float>(node.center.z - myPosZ);
            r2_thread = dxF * dxF + dyF * dyF + dzF * dzF + softSqF;
            const float sizeSq = node.size * node.size;
            // Always-direct в зоне softening — никогда не аппроксимируем близкие узлы.
            const bool nearField = (r2_thread < alwaysDirectR2);
            // Dehnen MAC: s² < θ²·(r - δ)². Если r ≤ δ, тело может лежать внутри
            // ball'а узла → принудительно открываем.
            const float r_node = sqrtf(r2_thread);
            const float effR   = r_node - node.delta;
            acceptMe = (!nearField) && (effR > 0.0f) && (sizeSq < thetaSq * effR * effR);
        }
        // Неактивный поток (idx >= numBodies) голосует "accept", чтобы не
        // заставлять варп бесполезно открывать узлы.
        const unsigned int openMask = __ballot_sync(0xffffffff, !acceptMe);

        if (openMask == 0u) {
            if (active) {
                const float invR  = rsqrtf(r2_thread);
                const float invR2 = invR * invR;
                const float invR3 = invR * invR2;
                const float scale = G_CONST_F * static_cast<float>(node.mass) * invR3;

                // Монополь: -GM/r³ · r̂  (в нашей конвенции d = source - body).
                double ax = static_cast<double>(scale * dxF);
                double ay = static_cast<double>(scale * dyF);
                double az = static_cast<double>(scale * dzF);

                // Квадрупольная поправка: a_q = (5G(d·Q·d)/2r⁷)·d - (G/r⁵)·Q·d
                const float Qxx = node.Qxx;
                const float Qyy = node.Qyy;
                const float Qxy = node.Qxy;
                const float Qxz = node.Qxz;
                const float Qyz = node.Qyz;
                const float Qzz = -(Qxx + Qyy);
                const float Qdx = Qxx*dxF + Qxy*dyF + Qxz*dzF;
                const float Qdy = Qxy*dxF + Qyy*dyF + Qyz*dzF;
                const float Qdz = Qxz*dxF + Qyz*dyF + Qzz*dzF;
                const float dQd = dxF*Qdx + dyF*Qdy + dzF*Qdz;
                const float invR5 = invR3 * invR2;
                const float invR7 = invR5 * invR2;
                const float c5 = G_CONST_F * invR5;
                const float c7 = 2.5f * G_CONST_F * invR7;
                ax += static_cast<double>(c7 * dQd * dxF - c5 * Qdx);
                ay += static_cast<double>(c7 * dQd * dyF - c5 * Qdy);
                az += static_cast<double>(c7 * dQd * dzF - c5 * Qdz);

                acc.x += ax;
                acc.y += ay;
                acc.z += az;
            }
            continue;
        }

        if (laneId == 0) {
            #pragma unroll
            for (int child = 0; child < 8; ++child) {
                const int next = node.children[child];
                if (next >= 0) {
                    const int top = sStackTop[warpId];
                    if (top < WARP_STACK_DEPTH) {
                        myStack[top] = next;
                        sStackTop[warpId] = top + 1;
                    } else {
                        atomicOr(statusFlags, STATUS_TRAVERSAL_STACK_OVERFLOW);
                    }
                }
            }
        }
        __syncwarp();
    }

    if (active) {
        outAcc[idx] = acc;
    }
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

void buildTreeAndComputeAcceleration(
    size_t numBodies,
    int blocks,
    float theta,
    double softeningSq,
    double3* outAcc) {
    resetBuildStateKernel<<<1, 1>>>(d_minMax, d_leafCount, d_nodeCounter);
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

    // Memsets не нужны: buildLeafMetadataKernel пишет leafStarts/leafEnds только
    // в позиции [0, leafCount), а initLeafNodesKernel читает их только в этом
    // же диапазоне. Stale-данные за пределами никогда не используются.
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

    // Quadrupole pass: сначала Q каждого листа от его тел, затем накопление
    // вверх по дереву через atomicAdd со Steiner-сдвигом.
    computeLeafQuadrupoleKernel<<<blocks, THREADS_PER_BLOCK>>>(
        d_nodes, d_posMass, d_sortedBodyIndices, d_leafNodeIndices, d_leafCount);
    CHECK_CUDA(cudaGetLastError());

    accumulateQuadrupoleKernel<<<blocks, THREADS_PER_BLOCK>>>(
        d_nodes, d_leafCodes, d_leafNodeIndices, d_leafCount);
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
    // statusFlags читаем не здесь, а один раз в simulationStep после всех substep'ов —
    // это убирает D→H sync stall на каждом substep (см. фикс #2).
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
    CHECK_CUDA(cudaMalloc(&d_maxAcc, sizeof(double)));
    CHECK_CUDA(cudaMallocHost(&h_maxAccPinned, sizeof(double)));
    CHECK_CUDA(cudaMallocHost(&h_statusFlagsPinned, sizeof(int)));

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
    if (d_maxAcc) cudaFree(d_maxAcc);
    if (h_maxAccPinned) cudaFreeHost(h_maxAccPinned);
    if (h_statusFlagsPinned) cudaFreeHost(h_statusFlagsPinned);

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
    d_maxAcc = nullptr;
    h_maxAccPinned = nullptr;
    h_statusFlagsPinned = nullptr;

    currentCapacity = 0;
    currentNodeCapacity = 0;
    deviceStateValid = false;
    accInitialized = false;
}

void simulationStepBarnesHutCUDA(std::vector<Object>& objs, double dt, bool pause, bool forceSync) {
    if (objs.empty()) {
        return;
    }
    // pause и forceSync — независимы: пауза значит "не интегрировать",
    // forceSync значит "синхронизировать host↔device". Не выходим из функции
    // на pause, иначе host не получит свежие данные с GPU когда renderer
    // запрашивает (state.pause + forceSync=true в realtime).

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

    if (dt > 0.0 && !pause) {
        // Один сброс status-флагов на весь outer-step. Между substep'ами они
        // атомарно ORятся; читаем после всего цикла.
        CHECK_CUDA(cudaMemset(d_statusFlags, 0, sizeof(int)));

        if (!accInitialized) {
            buildTreeAndComputeAcceleration(N, blocks, theta, softSq, d_acc);
            accInitialized = true;
        }

        // Adaptive substep с гарантированным завершением dt:
        // minStep = dt/cap → даже в худшем случае (плотный кластер,
        // η-формула просит h<<minStep) цикл завершится ровно за cap
        // итераций и пройдёт ровно dt — никакого drift'а симуляционного времени.
        const int    substepCap = std::max(1, physics::getSubstepCap());
        const double adaptiveEta = physics::getAdaptiveEta();
        const double minStep = dt / static_cast<double>(substepCap);
        double remaining = dt;
        int substeps = 0;
        while (remaining > 0.0 && substeps < substepCap) {
            // d_maxAcc сбрасываем через async-memset (без host stall), maxAccelKernel
            // ORит результат, читаем в pinned-память — это убирает overhead от
            // pageable-memory transfer.
            CHECK_CUDA(cudaMemsetAsync(d_maxAcc, 0, sizeof(double)));
            maxAccelKernel<<<blocks, THREADS_PER_BLOCK>>>(d_acc, static_cast<int>(N), d_maxAcc);
            CHECK_CUDA(cudaGetLastError());
            CHECK_CUDA(cudaMemcpy(h_maxAccPinned, d_maxAcc, sizeof(double), cudaMemcpyDeviceToHost));
            const double maxAcc = *h_maxAccPinned;

            double suggested = remaining;
            if (maxAcc > 0.0) {
                suggested = std::min(suggested, adaptiveEta / std::sqrt(maxAcc));
            }
            const double h = std::clamp(suggested, std::min(minStep, remaining), remaining);

            integratePositionKernel<<<blocks, THREADS_PER_BLOCK>>>(d_posMass, d_vel, d_acc, static_cast<int>(N), h);
            CHECK_CUDA(cudaGetLastError());

            buildTreeAndComputeAcceleration(N, blocks, theta, softSq, d_newAcc);
            updateVelocityKernel<<<blocks, THREADS_PER_BLOCK>>>(d_vel, d_acc, d_newAcc, static_cast<int>(N), h);
            CHECK_CUDA(cudaGetLastError());

            remaining -= h;
            ++substeps;
        }

        CHECK_CUDA(cudaMemcpy(h_statusFlagsPinned, d_statusFlags, sizeof(int), cudaMemcpyDeviceToHost));
        const int statusFlags = *h_statusFlagsPinned;
        if (statusFlags != 0) {
            if ((statusFlags & STATUS_BUILD_NODE_OVERFLOW) != 0) {
                throw std::runtime_error("Barnes-Hut CUDA tree capacity exhausted during build.");
            }
            if ((statusFlags & STATUS_TRAVERSAL_STACK_OVERFLOW) != 0) {
                throw std::runtime_error("Barnes-Hut CUDA traversal stack overflowed.");
            }
        }
    }

    // Хост-копию делаем только когда вызывающий явно просит — это либо кадр
    // на сохранение (headless), либо рендер (realtime), либо финальный cleanup.
    // На промежуточных headless-шагах без сохранения это 800+ KB D→H впустую.
    if (forceSync) {
        CHECK_CUDA(cudaMemcpy(h_posMass.data(), d_posMass, N * sizeof(double4), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(h_vel.data(), d_vel, N * sizeof(double3), cudaMemcpyDeviceToHost));
        CHECK_CUDA(cudaMemcpy(h_acc.data(), d_acc, N * sizeof(double3), cudaMemcpyDeviceToHost));
        for (size_t i = 0; i < N; ++i) {
            objs[i].position     = glm::dvec3(h_posMass[i].x, h_posMass[i].y, h_posMass[i].z);
            objs[i].velocity     = glm::dvec3(h_vel[i].x, h_vel[i].y, h_vel[i].z);
            objs[i].acceleration = glm::dvec3(h_acc[i].x, h_acc[i].y, h_acc[i].z);
        }
    }
}
