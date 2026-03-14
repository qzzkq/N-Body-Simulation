#include "brut_force.hpp"

#include <cmath>
#include <vector>

#include <object.hpp>
#include <physics.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

static constexpr double SOFTENING_AU = 1.0e6 * physics::METERS_TO_AU; // 1000 km

void computeAccelerations(const std::vector<glm::dvec3>& positions,
                          const std::vector<double>& masses,
                          std::vector<glm::dvec3>& accelerations) {
    const size_t n = positions.size();

    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n; ++i) {
        glm::dvec3 acc(0.0);
        const glm::dvec3 pi = positions[i];

        for (size_t j = 0; j < n; ++j) {
            if (i == j) {
                continue;
            }

            const glm::dvec3 delta = positions[j] - pi;
            const double distSq = delta.x * delta.x + delta.y * delta.y + delta.z * delta.z + SOFTENING_AU * SOFTENING_AU;

            const double invDist = 1.0 / std::sqrt(distSq);
            const double invDist3 = invDist * invDist * invDist;
            const double factor = physics::G * masses[j] * invDist3;

            acc += delta * factor;
        }

        accelerations[i] = acc;
    }
}

} // namespace

void simulationStepBrutForceCPU(std::vector<Object>& objs, double dt, bool pause) {
    if (pause || objs.empty()) {
        return;
    }

    const size_t n = objs.size();

    // Reused work buffers to minimize per-step allocations.
    static std::vector<double> masses;
    static std::vector<glm::dvec3> pos0;
    static std::vector<glm::dvec3> vel0;
    static std::vector<glm::dvec3> tmpPos;
    static std::vector<glm::dvec3> k1x, k2x, k3x, k4x;
    static std::vector<glm::dvec3> k1v, k2v, k3v, k4v;

    masses.resize(n);
    pos0.resize(n);
    vel0.resize(n);
    tmpPos.resize(n);
    k1x.resize(n);
    k2x.resize(n);
    k3x.resize(n);
    k4x.resize(n);
    k1v.resize(n);
    k2v.resize(n);
    k3v.resize(n);
    k4v.resize(n);

    for (size_t i = 0; i < n; ++i) {
        masses[i] = objs[i].mass;
        pos0[i] = objs[i].position;
        vel0[i] = objs[i].velocity;
    }

    // k1
    computeAccelerations(pos0, masses, k1v);
    for (size_t i = 0; i < n; ++i) {
        k1x[i] = vel0[i];
        tmpPos[i] = pos0[i] + 0.5 * dt * k1x[i];
    }

    // k2
    computeAccelerations(tmpPos, masses, k2v);
    for (size_t i = 0; i < n; ++i) {
        k2x[i] = vel0[i] + 0.5 * dt * k1v[i];
        tmpPos[i] = pos0[i] + 0.5 * dt * k2x[i];
    }

    // k3
    computeAccelerations(tmpPos, masses, k3v);
    for (size_t i = 0; i < n; ++i) {
        k3x[i] = vel0[i] + 0.5 * dt * k2v[i];
        tmpPos[i] = pos0[i] + dt * k3x[i];
    }

    // k4
    computeAccelerations(tmpPos, masses, k4v);
    for (size_t i = 0; i < n; ++i) {
        k4x[i] = vel0[i] + dt * k3v[i];
    }

    const double w = dt / 6.0;
    for (size_t i = 0; i < n; ++i) {
        objs[i].position = pos0[i] + w * (k1x[i] + 2.0 * k2x[i] + 2.0 * k3x[i] + k4x[i]);
        objs[i].velocity = vel0[i] + w * (k1v[i] + 2.0 * k2v[i] + 2.0 * k3v[i] + k4v[i]);
    }
}
