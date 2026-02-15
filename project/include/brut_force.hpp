#pragma once

#include <vector>
#include "object.hpp"

void simulationStepBrutForceCPU(std::vector<Object>& objs, float dt, bool pause, int iterations);
