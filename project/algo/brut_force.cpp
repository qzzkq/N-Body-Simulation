#include "brut_force.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include <object.hpp>
#include <physics.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

// O(N²) попарно с симметрией Ньютона. Параллелизуется через thread-local
// буфера; финальная свёртка — Neumaier-компенсированная.
void computeAccelerations(const std::vector<glm::dvec3>& positions,
                          const std::vector<double>& masses,
                          std::vector<glm::dvec3>& accelerations) {
    const size_t n = positions.size();
    accelerations.assign(n, glm::dvec3(0.0));
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
                if (rawDistSq <= 0.0 && soft2 <= 0.0) continue;
                const double invDist  = 1.0 / std::sqrt(distSq);
                const double invDist3 = invDist * invDist * invDist;
                const double common   = physics::G * invDist3;
                const glm::dvec3 ai = delta * (common * masses[j]);
                const glm::dvec3 aj = delta * (common * masses[static_cast<size_t>(i)]);
                acc[static_cast<size_t>(i)] += ai;
                acc[j]                       -= aj;
            }
        }
    }

    static std::vector<glm::dvec3> comp;
    if (comp.size() != n) comp.assign(n, glm::dvec3(0.0));
    else                  std::fill(comp.begin(), comp.end(), glm::dvec3(0.0));

    for (const auto& accThread : localAcc) {
        for (size_t i = 0; i < n; ++i) {
            const glm::dvec3 y = accThread[i] - comp[i];
            const glm::dvec3 t = accelerations[i] + y;
            comp[i]         = (t - accelerations[i]) - y;
            accelerations[i] = t;
        }
    }
#else
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            const glm::dvec3 delta = positions[j] - positions[i];
            const double rawDistSq = delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
            const double distSq = rawDistSq + soft2;
            if (rawDistSq <= 0.0 && soft2 <= 0.0) continue;
            const double invDist  = 1.0 / std::sqrt(distSq);
            const double invDist3 = invDist * invDist * invDist;
            const double common   = physics::G * invDist3;
            const glm::dvec3 ai = delta * (common * masses[j]);
            const glm::dvec3 aj = delta * (common * masses[i]);
            accelerations[i] += ai;
            accelerations[j] -= aj;
        }
    }
#endif
}

} // namespace

void simulationStepBrutForceCPU(std::vector<Object>& objs, double dt, bool pause, bool forceSync) {
    (void)forceSync;
    if (pause || objs.empty()) return;
    if (dt <= 0.0)              return;

    const size_t n = objs.size();

    static std::vector<double>     masses;
    static std::vector<glm::dvec3> x, v;
    static std::vector<glm::dvec3> aCur, aNext, tmpX;

    masses.resize(n);
    x.resize(n);
    v.resize(n);
    aCur.resize(n);
    aNext.resize(n);
    tmpX.resize(n);

    for (size_t i = 0; i < n; ++i) {
        masses[i] = objs[i].mass;
        x[i]      = objs[i].position;
        v[i]      = objs[i].velocity;
    }

    // Yoshida4 = три Velocity-Verlet с коэф. (w1, w0, w1); w0 < 0.
    // Фиксированный шаг — обязательное условие симплектичности интегратора.
    const double cbrt2 = std::cbrt(2.0);
    const double w1    = 1.0 / (2.0 - cbrt2);
    const double w0    = -cbrt2 * w1;
    const double yoshidaCoeffs[3] = { w1, w0, w1 };

    computeAccelerations(x, masses, aCur);

    auto verletSubstep = [&](double hs) {
        const double half_hs = 0.5 * hs;
        const double hs2     = hs * hs;
        for (size_t i = 0; i < n; ++i) {
            tmpX[i] = x[i] + v[i] * hs + 0.5 * aCur[i] * hs2;
        }
        computeAccelerations(tmpX, masses, aNext);
        for (size_t i = 0; i < n; ++i) {
            v[i]   += half_hs * (aCur[i] + aNext[i]);
            x[i]    = tmpX[i];
            aCur[i] = aNext[i]; // FSAL
        }
    };

    for (int s = 0; s < 3; ++s) {
        verletSubstep(yoshidaCoeffs[s] * dt);
    }

    for (size_t i = 0; i < n; ++i) {
        objs[i].position     = x[i];
        objs[i].velocity     = v[i];
        objs[i].acceleration = aCur[i];
    }
}
