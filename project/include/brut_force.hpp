/**
 * @file brut_force.hpp
 * @brief Алгоритм Брутфорс — прямой попарный расчёт гравитации O(N²).
 *
 * Для каждой пары тел (i, j) вычисляет силу Ньютона:
 * @f[
 *   \vec{F}_{ij} = G \frac{m_i m_j}{r_{ij}^2 + \varepsilon^2} \hat{r}_{ij}
 * @f]
 * где ε — параметр сглаживания (physics::getSofteningAU()).
 *
 * Интегратор: leapfrog (Verlet на скоростях).
 *
 * @par Производительность:
 *  - O(N²) взаимодействий за шаг.
 *  - Рекомендуется при N < 1 000.
 *  - При N = 10 000 на CPU: ~100 мс/шаг (зависит от железа).
 *
 * @par Точность:
 *  - Наиболее точный из трёх методов (нет апроксимации).
 *  - Погрешность определяется только шагом интегрирования dt.
 */

#pragma once
#include <vector>
#include "object.hpp"

void simulationStepBrutForceCPU(std::vector<Object>& objs, double dt, bool pause);
