#pragma once
#include <glm/glm.hpp>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include "object.hpp"

namespace physics {
    constexpr double PI = 3.14159265358979323846;

    // IAU Astronomical System of Units (normalized simulation units)
    constexpr double SOLAR_MASS_KG = 1.98847e30;
    constexpr double AU_METERS = 1.495978707e11;
    constexpr double YEAR_SECONDS = 3.15576e7;

    // SI -> normalized units conversion factors
    constexpr double MASS_TO_SOLAR = 1.0 / SOLAR_MASS_KG;
    constexpr double METERS_TO_AU = 1.0 / AU_METERS;
    constexpr double SECONDS_TO_YEARS = 1.0 / YEAR_SECONDS;
    constexpr double VELOCITY_TO_AU_PER_YEAR = METERS_TO_AU / SECONDS_TO_YEARS;

    // Gravitational constant in [AU^3 / (M_sun * year^2)]
    constexpr double G = 4.0 * PI * PI;
    // Must match integrator softening (currently 1000 km in AU)
    constexpr double SOFTENING_AU = 1.0e6 * METERS_TO_AU;

    constexpr double RADIUS_SCALE = 1.0;

    inline float calculateRadius(double mass_solar, double density_kg_m3) {
        if (density_kg_m3 <= 0.0) return 1.0f;
        const double mass_kg = mass_solar / MASS_TO_SOLAR;
        const double r_meters = std::cbrt((3.0 * mass_kg) / (4.0 * PI * density_kg_m3));
        return static_cast<float>((r_meters * METERS_TO_AU) / RADIUS_SCALE);
    }

    double calculateTotalEnergy(const std::vector<Object>& objs);

    void colorFromMass(std::vector<Object>& objs);
}
