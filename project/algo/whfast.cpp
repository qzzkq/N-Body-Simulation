#include "whfast.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

#include <object.hpp>
#include <physics.hpp>

namespace {

// Stumpff c_n(z): ряд для |z| ≤ 0.1, иначе половинная редукция и удвоение.
//   c_0(4z) = 2c_0² - 1,  c_1(4z) = c_0·c_1,
//   c_2(4z) = c_1²/2,     c_3(4z) = (c_2 + c_0·c_3)/4
void stumpff(double z, double& c0, double& c1, double& c2, double& c3) {
    int n = 0;
    while (std::abs(z) > 0.1 && n < 32) {
        z *= 0.25;
        ++n;
    }
    c3 = (1.0 / 6.0)
       * (1.0 - z * (1.0 / 20.0)
                * (1.0 - z * (1.0 / 42.0)
                         * (1.0 - z * (1.0 / 72.0)
                                  * (1.0 - z / 110.0))));
    c2 = 0.5
       * (1.0 - z * (1.0 / 12.0)
                * (1.0 - z * (1.0 / 30.0)
                         * (1.0 - z * (1.0 / 56.0)
                                  * (1.0 - z / 90.0))));
    c1 = 1.0 - z * c3;
    c0 = 1.0 - z * c2;
    while (n > 0) {
        c3 = (c2 + c0 * c3) * 0.25;
        c2 = (c1 * c1) * 0.5;
        c1 = c0 * c1;
        c0 = 2.0 * c0 * c0 - 1.0;
        --n;
    }
}

// Newton по χ для уравнения Баттина:
//   √μ·dt = r0·χ·c1 + (σ0/√μ)·χ²·c2 + χ³·c3,   F'(χ) = r(χ).
void kepler_step(glm::dvec3& r, glm::dvec3& v, double mu, double dt) {
    if (dt == 0.0) return;
    const double r0 = glm::length(r);
    if (r0 == 0.0) return;
    const double v0_sq  = glm::dot(v, v);
    const double sigma0 = glm::dot(r, v);
    const double alpha  = 2.0 / r0 - v0_sq / mu;
    const double smu    = std::sqrt(mu);
    const double smu_dt = smu * dt;

    double X = smu * dt / r0;
    double c0 = 0.0, c1 = 0.0, c2 = 0.0, c3 = 0.0;

    for (int iter = 0; iter < 32; ++iter) {
        const double z = alpha * X * X;
        stumpff(z, c0, c1, c2, c3);
        const double X_c1  = X * c1;
        const double X2_c2 = X * X * c2;
        const double X3_c3 = X * X * X * c3;
        const double F  = r0 * X_c1 + (sigma0 / smu) * X2_c2 + X3_c3 - smu_dt;
        const double Fp = r0 * c0   + (sigma0 / smu) * X_c1  + X2_c2;
        const double dX = F / Fp;
        X -= dX;
        if (std::abs(dX) < 1e-14 * std::max(1.0, std::abs(X))) break;
    }

    {
        const double z = alpha * X * X;
        stumpff(z, c0, c1, c2, c3);
    }

    const double f = 1.0 - (X * X / r0) * c2;
    const double g = dt - (X * X * X / smu) * c3;

    const glm::dvec3 r_new = f * r + g * v;
    const double r_new_len = glm::length(r_new);
    if (r_new_len == 0.0) return;

    const double f_dot = -smu * X * c1 / (r_new_len * r0);
    const double g_dot = 1.0 - (X * X / r_new_len) * c2;

    const glm::dvec3 v_new = f_dot * r + g_dot * v;
    r = r_new;
    v = v_new;
}

void interaction_kick(std::vector<glm::dvec3>& Q,
                      std::vector<glm::dvec3>& V,
                      const std::vector<double>& m,
                      size_t central, double h) {
    const size_t n = Q.size();
    for (size_t i = 0; i < n; ++i) {
        if (i == central) continue;
        glm::dvec3 ai(0.0);
        for (size_t j = 0; j < n; ++j) {
            if (j == central || j == i) continue;
            const glm::dvec3 d = Q[j] - Q[i];
            const double r2 = glm::dot(d, d);
            if (r2 == 0.0) continue;
            const double inv_r  = 1.0 / std::sqrt(r2);
            const double inv_r3 = inv_r * inv_r * inv_r;
            ai += (physics::G * m[j] * inv_r3) * d;
        }
        V[i] += h * ai;
    }
}

void sun_jump(std::vector<glm::dvec3>& Q,
              const std::vector<glm::dvec3>& V,
              const std::vector<double>& m,
              size_t central, double h) {
    const size_t n = Q.size();
    glm::dvec3 P_planets(0.0);
    for (size_t i = 0; i < n; ++i) {
        if (i == central) continue;
        P_planets += m[i] * V[i];
    }
    const glm::dvec3 shift = (h / m[central]) * P_planets;
    for (size_t i = 0; i < n; ++i) {
        if (i == central) continue;
        Q[i] += shift;
    }
}

} // namespace

void simulationStepWHFastCPU(std::vector<Object>& objs, double dt, bool pause, bool forceSync) {
    (void)forceSync;
    if (pause || objs.empty()) return;
    if (dt == 0.0)              return;
    if (objs.size() < 2)        return;

    const size_t n = objs.size();

    size_t central = 0;
    for (size_t i = 1; i < n; ++i)
        if (objs[i].mass > objs[central].mass) central = i;

    static std::vector<double>     m;
    static std::vector<glm::dvec3> Q, V;
    m.resize(n);
    Q.resize(n);
    V.resize(n);

    glm::dvec3 r_com(0.0), v_com(0.0);
    double M_tot = 0.0;
    for (size_t i = 0; i < n; ++i) {
        m[i]   = objs[i].mass;
        M_tot += m[i];
        r_com += m[i] * objs[i].position;
        v_com += m[i] * objs[i].velocity;
    }
    if (M_tot <= 0.0 || m[central] <= 0.0) return;
    r_com /= M_tot;
    v_com /= M_tot;

    const glm::dvec3 r_central_in = objs[central].position;
    for (size_t i = 0; i < n; ++i) {
        Q[i] = objs[i].position - r_central_in;
        V[i] = objs[i].velocity - v_com;
    }

    const double mu = physics::G * m[central];
    const double h  = dt;

    // WHFast 2-го порядка: I(h/2) ∘ J(h/2) ∘ K(h) ∘ J(h/2) ∘ I(h/2).
    interaction_kick(Q, V, m, central, 0.5 * h);
    sun_jump        (Q, V, m, central, 0.5 * h);
    for (size_t i = 0; i < n; ++i) {
        if (i == central) continue;
        kepler_step(Q[i], V[i], mu, h);
    }
    sun_jump        (Q, V, m, central, 0.5 * h);
    interaction_kick(Q, V, m, central, 0.5 * h);

    // Восстановление позиции/скорости центрального тела из ∑ m_j·r_j = 0.
    glm::dvec3 r_central_bary(0.0);
    glm::dvec3 V_central     (0.0);
    for (size_t i = 0; i < n; ++i) {
        if (i == central) continue;
        r_central_bary -= m[i] * Q[i];
        V_central      -= m[i] * V[i];
    }
    r_central_bary /= M_tot;
    V_central      /= m[central];

    const glm::dvec3 r_com_new = r_com + v_com * dt;

    for (size_t i = 0; i < n; ++i) {
        glm::dvec3 r_bary, v_bary;
        if (i == central) {
            r_bary = r_central_bary;
            v_bary = V_central;
        } else {
            r_bary = r_central_bary + Q[i];
            v_bary = V[i];
        }
        objs[i].position = r_com_new + r_bary;
        objs[i].velocity = v_com    + v_bary;
    }

    // Полное Newton-поле ускорений — только для визуализации.
    for (size_t i = 0; i < n; ++i) {
        glm::dvec3 a(0.0);
        for (size_t j = 0; j < n; ++j) {
            if (j == i) continue;
            const glm::dvec3 d = objs[j].position - objs[i].position;
            const double r2 = glm::dot(d, d);
            if (r2 == 0.0) continue;
            const double inv_r  = 1.0 / std::sqrt(r2);
            const double inv_r3 = inv_r * inv_r * inv_r;
            a += (physics::G * m[j] * inv_r3) * d;
        }
        objs[i].acceleration = a;
    }
}
