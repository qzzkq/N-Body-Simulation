#pragma once
#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "object.hpp"

#define CHECK_CUDA(call)                                                        \
    do {                                                                        \
        cudaError_t err__ = (call);                                             \
        if (err__ != cudaSuccess) {                                             \
            std::fprintf(stderr, "CUDA error at %s:%d in %s: %s\n",             \
                         __FILE__, __LINE__, #call, cudaGetErrorString(err__)); \
            std::abort();                                                       \
        }                                                                       \
    } while (0)

void InitBarnesHutCUDA(size_t maxObjects);

void CleanupBarnesHutCUDA();

void simulationStepBarnesHutCUDA(std::vector<Object>& objs, double dt, bool pause, bool forceSync = false);
