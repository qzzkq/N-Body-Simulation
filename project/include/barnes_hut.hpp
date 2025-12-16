// barnes_hut_step.hpp
#pragma once

#include <vector>
#include "object.hpp"

void simulationStepBarnesHutCPU(std::vector<Object>& objs, float dt, bool pause, int iterations);
