/**
 * @file graphic_state.hpp
 * @brief Визуальное состояние тела: цвет и история позиций (след).
 *
 * GraphicState специально отделён от Object, чтобы физический слой
 * не зависел от слоя визуализации. Массив GraphicState живёт рядом
 * с массивом Object и индексируется синхронно с ним.
 */

<<<<<<< gui
#include <array>
#include <glm/glm.hpp>

static constexpr std::size_t MAX_TRAIL_LENGTH = 150;

=======
#pragma once
#include <deque>
#include <glm/glm.hpp>

/**
 * @struct GraphicState
 * @brief Визуальные данные одного тела: цвет и «хвост» траектории.
 *
 * Обновляется на каждом кадре рендеринга и при воспроизведении
 * записанной симуляции (Replay). Хранится отдельно от Object,
 * чтобы сохранить чистоту физического слоя.
 *
 * @note В режиме Replay хвост строится заново при перемотке назад —
 *       метод resetTrail() очищает @c trail.
 */
>>>>>>> main
struct GraphicState {

    /**
     * @brief RGBA-цвет тела [0, 1].
     *
     * Устанавливается при загрузке из файла (если .h5 содержит секцию
     * @c body_colors) или вычисляется через physics::colorFromMass().
     * По умолчанию — белый (1, 1, 1, 1).
     */
    glm::vec4 color = glm::vec4(1.0f);

<<<<<<< gui
    std::array<glm::vec3, MAX_TRAIL_LENGTH> trail{};
    int trailHead = 0;  // ring buffer write pointer
    int trailSize = 0;  // [0, MAX_TRAIL_LENGTH]

    // Логический индекс: 0 = самая старая точка, trailSize-1 = самая новая
    const glm::vec3& trailAt(int i) const {
        int start = (trailHead - trailSize + static_cast<int>(MAX_TRAIL_LENGTH)) % static_cast<int>(MAX_TRAIL_LENGTH);
        return trail[static_cast<std::size_t>((start + i) % static_cast<int>(MAX_TRAIL_LENGTH))];
    }
=======
    /**
     * @brief Кольцевой буфер последних позиций тела — «хвост» траектории.
     *
     * Используется Renderer::drawTrails() для рисования GL_LINE_STRIP.
     * Ёмкость ограничена MAX_TRAIL_LENGTH: при переполнении самая старая
     * точка удаляется с переднего конца.
     *
     * Точки добавляются только если тело переместилось дальше чем
     * @c distanceThreshold (0.2 AU), чтобы не засорять буфер.
     */
    std::deque<glm::vec3> trail;
>>>>>>> main

    /**
     * @brief Добавляет точку в хвост, если тело переместилось достаточно далеко.
     *
     * Алгоритм:
     *  1. Если @p force == @c true или буфер пуст — добавляет всегда.
     *  2. Иначе добавляет только если расстояние от последней точки
     *     превышает @c distanceThreshold (0.2 AU).
     *  3. Если буфер достиг MAX_TRAIL_LENGTH — удаляет самую старую точку.
     *
     * Вызывается из главного цикла симуляции (@c main.cpp) и Replay после
     * каждого шага интеграции или кадра воспроизведения.
     *
     * @param position  Текущая позиция тела [AU], double → приводится к float.
     * @param force     @c true — добавить точку независимо от расстояния.
     */
    void updateTrail(const glm::dvec3& position, bool force = false) {
<<<<<<< gui
        const glm::vec3 pos32(position);

        if (!force && trailSize > 0) {
            constexpr float distanceThreshold = 0.2f;
            const glm::vec3& last = trail[static_cast<std::size_t>(
                (trailHead - 1 + static_cast<int>(MAX_TRAIL_LENGTH)) % static_cast<int>(MAX_TRAIL_LENGTH))];
            if (glm::distance(pos32, last) <= distanceThreshold) return;
        }

        trail[static_cast<std::size_t>(trailHead)] = pos32;
        trailHead = (trailHead + 1) % static_cast<int>(MAX_TRAIL_LENGTH);
        if (trailSize < static_cast<int>(MAX_TRAIL_LENGTH)) ++trailSize;
=======
        if (force || trail.empty()) {
            trail.push_back(glm::vec3(position));
            return;
        }
        constexpr float distanceThreshold = 0.2f;
        if (glm::distance(glm::vec3(position), trail.back()) > distanceThreshold) {
            trail.push_back(glm::vec3(position));
        }
        if (trail.size() > MAX_TRAIL_LENGTH) {
            trail.pop_front();
        }
>>>>>>> main
    }

    /**
     * @brief Полностью очищает хвост траектории.
     *
     * Вызывается при перемотке назад в Replay (обнаружение
     * уменьшения frameA), чтобы хвост не рисовался «задом наперёд».
     */
    void resetTrail() {
        trailHead = 0;
        trailSize = 0;
    }
<<<<<<< gui
=======

    /**
     * @brief Максимальная длина хвоста траектории (количество точек).
     *
     * При 60 FPS и типичном шаге 1 AU / раз в 3 кадра
     * хвост отображает примерно 7–8 лет орбиты.
     */
    static constexpr std::size_t MAX_TRAIL_LENGTH = 150;
>>>>>>> main
};
