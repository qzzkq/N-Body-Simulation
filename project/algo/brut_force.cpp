#include "brut_force.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

#include <object.hpp>
#include <physics.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

struct AccelStats {
    double maxAcc = 0.0;
    double minDist = std::numeric_limits<double>::max();
};

void computeAccelerations(const std::vector<glm::dvec3>& positions,
                          const std::vector<double>& masses,
                          std::vector<glm::dvec3>& accelerations,
                          AccelStats* stats = nullptr) {
    const size_t n = positions.size();
    accelerations.assign(n, glm::dvec3(0.0));
    double minDistLocal = std::numeric_limits<double>::max();
    const double soft = physics::getSofteningAU();
    const double soft2 = soft * soft;

#ifdef _OPENMP
    const int threadCount = std::max(1, omp_get_max_threads());
    static std::vector<std::vector<glm::dvec3>> localAcc;
    if (localAcc.size() != static_cast<std::size_t>(threadCount) ||
        (!localAcc.empty() && localAcc.front().size() != n)) {
        localAcc.assign(
            static_cast<std::size_t>(threadCount),
            std::vector<glm::dvec3>(n, glm::dvec3(0.0))
        );
    } else {
        for (auto& accArr : localAcc) {
            std::fill(accArr.begin(), accArr.end(), glm::dvec3(0.0));
        }
    }

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        auto& acc = localAcc[static_cast<std::size_t>(tid)];

        #pragma omp for schedule(static)
        for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(n); ++i) {
            #pragma omp simd
            for (size_t j = static_cast<size_t>(i) + 1; j < n; ++j) {
                const glm::dvec3 delta = positions[j] - positions[static_cast<size_t>(i)];
                const double rawDistSq = delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
                const double distSq = rawDistSq + soft2;
                if (rawDistSq > 0.0) {
                    const double rawDist = std::sqrt(rawDistSq);
                    if (rawDist < minDistLocal) minDistLocal = rawDist;
                } else if (soft2 <= 0.0) {
                    continue;
                }
                const double invDist = 1.0 / std::sqrt(distSq);
                const double invDist3 = invDist * invDist * invDist;

                const double common = physics::G * invDist3;
                const glm::dvec3 ai = delta * (common * masses[j]);
                const glm::dvec3 aj = delta * (common * masses[static_cast<size_t>(i)]);

                acc[static_cast<size_t>(i)] += ai;
                acc[j] -= aj;
            }
        }
    }

    for (const auto& accThread : localAcc) {
        for (size_t i = 0; i < n; ++i) {
            accelerations[i] += accThread[i];
        }
    }
#else
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            const glm::dvec3 delta = positions[j] - positions[i];
            const double rawDistSq = delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
            const double distSq = rawDistSq + soft2;
            if (rawDistSq > 0.0) {
                const double rawDist = std::sqrt(rawDistSq);
                if (rawDist < minDistLocal) minDistLocal = rawDist;
            } else if (soft2 <= 0.0) {
                continue;
            }
            const double invDist = 1.0 / std::sqrt(distSq);
            const double invDist3 = invDist * invDist * invDist;

            const double common = physics::G * invDist3;
            const glm::dvec3 ai = delta * (common * masses[j]);
            const glm::dvec3 aj = delta * (common * masses[i]);

            accelerations[i] += ai;
            accelerations[j] -= aj;
        }
    }
#endif

    if (stats != nullptr) {
        double maxA = 0.0;
        for (const auto& a : accelerations) {
            maxA = std::max(maxA, glm::length(a));
        }
        stats->maxAcc = maxA;
        stats->minDist = minDistLocal;
    }
}

} // namespace

void simulationStepBrutForceCPU(std::vector<Object>& objs, double dt, bool pause) {
    if (pause || objs.empty()) {
        return;
    }
    if (dt <= 0.0) {
        return;
    }

    const size_t n = objs.size();

    // Reused work buffers to minimize per-step allocations.
    static std::vector<double> masses;
    static std::vector<glm::dvec3> x, v;
    static std::vector<glm::dvec3> k1x, k2x, k3x, k4x;
    static std::vector<glm::dvec3> k1v, k2v, k3v, k4v;
    static std::vector<glm::dvec3> tmpX, tmpV;

    masses.resize(n);
    x.resize(n);
    v.resize(n);
    k1x.resize(n); k2x.resize(n); k3x.resize(n); k4x.resize(n);
    k1v.resize(n); k2v.resize(n); k3v.resize(n); k4v.resize(n);
    tmpX.resize(n);
    tmpV.resize(n);

    for (size_t i = 0; i < n; ++i) {
        masses[i] = objs[i].mass;
        x[i] = objs[i].position;
        v[i] = objs[i].velocity;
    }

    double remaining = dt;
    constexpr int MAX_SUBSTEPS = 1024;
    int substeps = 0;
    while (remaining > 0.0 && substeps < MAX_SUBSTEPS) {
        AccelStats st;
        computeAccelerations(x, masses, k1v, &st);
        for (size_t i = 0; i < n; ++i) {
            k1x[i] = v[i];
        }

        double suggested = remaining;
        if (st.maxAcc > 0.0) {
            suggested = std::min(suggested, 0.02 / std::sqrt(st.maxAcc));
        }
        if (st.minDist < std::numeric_limits<double>::max()) {
            const double maxMass = *std::max_element(masses.begin(), masses.end());
            const double dyn = 0.1 * std::sqrt((st.minDist * st.minDist * st.minDist) / std::max(physics::G * maxMass, 1e-30));
            suggested = std::min(suggested, dyn);
        }
        const double minStep = dt / 4096.0;
        double h = std::clamp(suggested, minStep, remaining);

        for (size_t i = 0; i < n; ++i) {
            tmpX[i] = x[i] + 0.5 * h * k1x[i];
            tmpV[i] = v[i] + 0.5 * h * k1v[i];
        }
        computeAccelerations(tmpX, masses, k2v, nullptr);
        for (size_t i = 0; i < n; ++i) {
            k2x[i] = tmpV[i];
        }

        for (size_t i = 0; i < n; ++i) {
            tmpX[i] = x[i] + 0.5 * h * k2x[i];
            tmpV[i] = v[i] + 0.5 * h * k2v[i];
        }
        computeAccelerations(tmpX, masses, k3v, nullptr);
        for (size_t i = 0; i < n; ++i) {
            k3x[i] = tmpV[i];
        }

        for (size_t i = 0; i < n; ++i) {
            tmpX[i] = x[i] + h * k3x[i];
            tmpV[i] = v[i] + h * k3v[i];
        }
        computeAccelerations(tmpX, masses, k4v, nullptr);
        for (size_t i = 0; i < n; ++i) {
            k4x[i] = tmpV[i];
        }

        const double w = h / 6.0;
        for (size_t i = 0; i < n; ++i) {
            x[i] += w * (k1x[i] + 2.0 * k2x[i] + 2.0 * k3x[i] + k4x[i]);
            v[i] += w * (k1v[i] + 2.0 * k2v[i] + 2.0 * k3v[i] + k4v[i]);
        }

        remaining -= h;
        ++substeps;
    }

    // Save state back to objects
    computeAccelerations(x, masses, k1v, nullptr);
    for (size_t i = 0; i < n; ++i) {
        objs[i].position = x[i];
        objs[i].velocity = v[i];
        objs[i].acceleration = k1v[i];
    }
}
