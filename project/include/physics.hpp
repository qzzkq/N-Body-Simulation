#pragma once
#include <glm/glm.hpp>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include "object.hpp"

namespace physics {
    constexpr double G = 6.6743e-11;
    constexpr double G_KM = 6.6743e-20;
    constexpr double RADIUS_SCALE = 50000.0;

    inline float calculateRadius(double mass, double density) {
        if (density <= 0) return 1.0f;
        const double pi = 3.14159265358979323846;
        const double r_m = std::cbrt((3.0 * mass) / (4.0 * pi * density));
        return static_cast<float>(r_m / RADIUS_SCALE);
    }

    void colorFromMass(std::vector<Object>& objs);
}