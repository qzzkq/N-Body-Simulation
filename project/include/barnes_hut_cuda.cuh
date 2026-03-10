#pragma once
#include <vector>
#include "object.hpp"

void InitBarnesHutCUDA(size_t maxObjects);

void CleanupBarnesHutCUDA();

void simulationStepBarnesHutCUDA(std::vector<Object>& objs, float dt, bool pause);