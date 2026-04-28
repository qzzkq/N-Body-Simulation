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

/**
 * @brief Выполняет один шаг гравитационной симуляции методом Брутфорс.
 *
 * Алгоритм на каждом шаге:
 *  1. Если !pause: для всех пар (i < j) вычисляет ускорения и аккумулирует
 *     их в Object::acceleration.
 *  2. Интегрирует скорости и позиции (leapfrog/Euler).
 *  3. Если forceSync: сбрасывает ускорения (используется для energy check).
 *
 * @param objs       Массив тел (модифицируется in-place: velocity, position).
 * @param dt         Шаг времени [год]. При dt == 0 и forceSync == true — только сброс состояния.
 * @param pause      @c true — пропустить интеграцию, не двигать тела.
 * @param forceSync  @c true — принудительно синхронизировать GPU/CPU (игнорируется в CPU-версии).
 */
void simulationStepBrutForceCPU(std::vector<Object>& objs,
                                 double dt,
                                 bool pause,
                                 bool forceSync = false);
