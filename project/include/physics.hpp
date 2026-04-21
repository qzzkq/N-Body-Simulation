/**
 * @file physics.hpp
 * @brief Физические константы, единицы измерения и вспомогательные функции симуляции.
 *
 * Все расчёты ведутся в нормализованной системе единиц IAU:
 * | Величина | Единица          |
 * |----------|------------------|
 * | Длина    | AU               |
 * | Масса    | M☉               |
 * | Время    | юлианский год    |
 * | G        | 4π² AU³/(M☉·год²)|
 *
 * Помимо констант модуль предоставляет:
 *  - Runtime-настраиваемые параметры интегратора (softening, theta).
 *  - Расчёт полной механической энергии системы (для верификации).
 *  - Окраску тел по массе для визуализации.
 */

#pragma once
#include <glm/glm.hpp>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include "object.hpp"
#include "graphic_state.hpp"

/**
 * @namespace physics
 * @brief Физические константы и утилиты для N-body симуляции.
 */
namespace physics {

    constexpr double PI = 3.14159265358979323846; ///< Математическая константа π.

    // ── Физические константы (SI) ───────────────────────────────────────────
    constexpr double SOLAR_MASS_KG  = 1.98847e30;       ///< Масса Солнца [кг].
    constexpr double AU_METERS      = 1.495978707e11;    ///< 1 AU в метрах.
    constexpr double YEAR_SECONDS   = 3.15576e7;         ///< 1 юлианский год в секундах.

    // ── Коэффициенты перевода SI → нормализованные единицы ─────────────────
    constexpr double MASS_TO_SOLAR            = 1.0 / SOLAR_MASS_KG;  ///< кг → M☉.
    constexpr double METERS_TO_AU             = 1.0 / AU_METERS;       ///< м → AU.
    constexpr double SECONDS_TO_YEARS         = 1.0 / YEAR_SECONDS;    ///< с → год.
    constexpr double VELOCITY_TO_AU_PER_YEAR  = METERS_TO_AU / SECONDS_TO_YEARS; ///< (м/с) → (AU/год).

    /**
     * @brief Гравитационная постоянная в нормализованных единицах.
     *
     * G = 4π²  [AU³ / (M☉ · год²)].
     * Выбрана так, что период обращения Земли вокруг Солнца на расстоянии
     * 1 AU при массе 1 M☉ равен ровно 1 году.
     */
    constexpr double G = 4.0 * PI * PI;

    constexpr double DEFAULT_SOFTENING_AU       = 0.0; ///< Параметр сглаживания гравитации по умолчанию [AU]. 0 = жёсткая физика.
    constexpr double DEFAULT_BARNES_HUT_THETA   = 0.5; ///< Параметр точности Barnes-Hut θ по умолчанию. Меньше → точнее, но медленнее.
    constexpr double RADIUS_SCALE               = 1.0; ///< Масштабный коэффициент радиуса при расчёте из плотности.

    /**
     * @brief Вычисляет физический радиус тела из его массы и плотности.
     *
     * Тело считается однородным шаром:
     * @f$ r = \left(\frac{3m}{4\pi\rho}\right)^{1/3} @f$
     *
     * Результат переводится в AU и масштабируется коэффициентом RADIUS_SCALE.
     *
     * @param mass_solar      Масса тела [M☉].
     * @param density_kg_m3   Средняя плотность [кг/м³]. При density ≤ 0 возвращает 1.0f.
     * @return Визуальный радиус [AU].
     */
    inline float calculateRadius(double mass_solar, double density_kg_m3) {
        if (density_kg_m3 <= 0.0) return 1.0f;
        const double mass_kg   = mass_solar / MASS_TO_SOLAR;
        const double r_meters  = std::cbrt((3.0 * mass_kg) / (4.0 * PI * density_kg_m3));
        return static_cast<float>((r_meters * METERS_TO_AU) / RADIUS_SCALE);
    }

    /**
     * @brief Вычисляет полную механическую энергию системы.
     *
     * E = Eкин + Eпот = Σ ½mᵢv²ᵢ  −  Σ G·mᵢ·mⱼ / rᵢⱼ
     *
     * Используется для верификации численной схемы: при идеальном интеграторе
     * энергия должна сохраняться. Сложность O(N²) — вызывать редко.
     *
     * @param objs Массив всех тел системы.
     * @return Полная механическая энергия [M☉·AU²/год²].
     */
    double calculateTotalEnergy(const std::vector<Object>& objs);

    /**
     * @brief Возвращает текущий параметр сглаживания гравитации.
     *
     * Softening ε предотвращает сингулярность при близком сближении:
     * F ~ G·m₁·m₂ / (r² + ε²).
     *
     * @return Значение ε [AU].
     */
    double getSofteningAU();

    /**
     * @brief Устанавливает параметр сглаживания гравитации.
     * @param softeningAu Новое значение ε [AU]. 0 — отключить сглаживание.
     */
    void setSofteningAU(double softeningAu);

    /**
     * @brief Возвращает текущий параметр точности Barnes-Hut θ.
     *
     * θ определяет критерий «далеко/близко» при обходе октодерева:
     *  s / d < θ, где s — размер ячейки, d — расстояние до тела.
     *  - θ = 0: эквивалентно Brute Force (все узлы раскрываются).
     *  - θ = 0.5: баланс скорости и точности (по умолчанию).
     *  - θ = 1.0+: быстро, но заметная погрешность.
     *
     * @return Текущее значение θ.
     */
    double getBarnesHutTheta();

    /**
     * @brief Устанавливает параметр точности Barnes-Hut.
     * @param theta Новое значение θ. Рекомендуемый диапазон: [0.3, 1.0].
     */
    void setBarnesHutTheta(double theta);

    /**
     * @brief Назначает цвета телам по логарифму их массы (устаревший API).
     *
     * @deprecated Используйте перегрузку с @c GraphicState.
     *             Цвет записывается как @c glm::vec4 напрямую в объект,
     *             что нарушает разделение физики и визуализации.
     * @param objs Массив тел (модифицируется).
     */
    void colorFromMass(std::vector<Object>& objs);

    /**
     * @brief Назначает цвета телам по логарифму их массы.
     *
     * Цветовая шкала (лог-нормирование):
     *  - малые массы (планеты): голубой/белый
     *  - средние (звёзды): жёлтый/оранжевый
     *  - гигантские (центральные): красный/ярко-белый
     *
     * @param objs     Массив тел (только для чтения масс).
     * @param graphics Массив визуальных состояний, куда записываются цвета.
     *                 Размер должен совпадать с @p objs.
     */
    void colorFromMass(const std::vector<Object>& objs, std::vector<GraphicState>& graphics);

} // namespace physics
