// barnes_hut_step.hpp
#pragma once

#include <vector>
#include "object.hpp"

void simulationStep(std::vector<Object>& objs, float dt, bool pause);
