/**
 * @file object.hpp
 * @brief Определение класса Object — базовой физической сущности симуляции.
 *
 * Каждый объект представляет собой тело с массой, позицией, скоростью
 * и ускорением в нормализованных астрономических единицах (AU, M☉, лет).
 *
 * Система единиц:
 *  - Длина:  астрономические единицы (AU)
 *  - Масса:  солнечные массы (M☉)
 *  - Время:  юлианские годы (365.25 суток)
 *  - G = 4π²  [AU³ / (M☉ · год²)]
 */

#pragma once
#include <glm/glm.hpp>
#include <optional>
#include <string>

/**
 * @class Object
 * @brief Физическое тело N-body симуляции.
 *
 * Хранит кинематические и физические характеристики тела.
 * Ускорение обнуляется и пересчитывается на каждом шаге интеграции
 * выбранным алгоритмом (BruteForce / Barnes-Hut CPU / Barnes-Hut CUDA).
 *
 * @par Пример создания:
 * @code
 * Object sun(
 *     glm::dvec3(0.0),        // позиция в центре масс
 *     glm::dvec3(0.0),        // скорость покоя
 *     1.0,                    // масса: 1 M☉
 *     1408.0,                 // плотность, кг/м³ (Солнце ≈ 1408)
 *     std::nullopt            // радиус вычислится из массы и плотности
 * );
 * @endcode
 */
class Object {
public:
    glm::dvec3 position;                    ///< Позиция в пространстве [AU].
    glm::dvec3 velocity;                    ///< Скорость [AU / год].
    glm::dvec3 acceleration = glm::dvec3(0.0); ///< Гравитационное ускорение [AU / год²]. Обнуляется на каждом шаге.
    double mass;                            ///< Масса тела [M☉].
    double density;                         ///< Плотность [кг / м³], используется для расчёта радиуса.

    float  radius;                          ///< Визуальный радиус тела [AU]. Вычисляется из mass и density.
    std::string name;                       ///< Имя тела (для именованных сценариев, например "Sun", "Earth").

    bool Initalizing = false; ///< @c true пока тело ещё «создаётся» интерактивно (клик мышью, не завершён).
    bool Launched    = false; ///< @c true после того как тело было запущено с заданной скоростью.
    bool target      = false; ///< @c true если тело выбрано как цель камеры.

    /**
     * @brief Конструктор тела.
     *
     * @param initPosition  Начальная позиция [AU].
     * @param initVelocity  Начальная скорость [AU/год].
     * @param massOpt       Масса [M☉]. Если @c nullopt — устанавливается значение по умолчанию.
     * @param densityOpt    Плотность [кг/м³]. Если @c nullopt — используется плотность воды (1000).
     * @param radiusOpt     Явный радиус [AU]. Если @c nullopt — вычисляется через physics::calculateRadius().
     */
    Object(glm::dvec3 initPosition,
           glm::dvec3 initVelocity,
           std::optional<double> massOpt,
           std::optional<double> densityOpt,
           std::optional<double> radiusOpt);

    /**
     * @brief Интегрирует позицию методом Эйлера.
     *
     * Обновляет @c position += velocity * deltaTime.
     * Вызывается после того как алгоритм (BruteForce/Barnes-Hut)
     * уже заполнил @c acceleration.
     *
     * @param deltaTime Шаг времени [год].
     */
    void UpdatePos(double deltaTime);

    /**
     * @brief Накапливает ускорение от одного взаимодействия.
     *
     * Прибавляет (ax, ay, az) * dt к @c velocity.
     * Вызывается внутри алгоритмов интеграции на каждое парное взаимодействие.
     *
     * @param x          Компонента ускорения по X [AU/год²].
     * @param y          Компонента ускорения по Y [AU/год²].
     * @param z          Компонента ускорения по Z [AU/год²].
     * @param deltaTime  Шаг времени [год].
     */
    void accelerate(double x, double y, double z, double deltaTime);

    /**
     * @brief Возвращает текущую позицию тела.
     * @return Вектор позиции [AU].
     */
    glm::dvec3 GetPos() const;
};
