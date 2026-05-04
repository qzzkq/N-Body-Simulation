/**
 * @file camera.hpp
 * @brief Камера наблюдателя в 3D-пространстве симуляции.
 *
 * Реализует модель «свободной камеры» (free-look / fly-cam):
 * движение по шести направлениям (WASD + QE) и вращение мышью
 * через углы Эйлера yaw / pitch.
 *
 * Матрица вида строится функцией glm::lookAt() на основе
 * текущих pos, front, up.
 */

#ifndef CAMERAINC
#define CAMERAINC

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

/**
 * @struct Camera
 * @brief Состояние камеры наблюдателя.
 *
 * Обновляется в Control::updateCameraFromKeys() каждый кадр.
 * Матрица вида вычисляется лениво через getViewMatrix() и передаётся
 * в Renderer::renderFrame() для построения MVP-матриц.
 *
 * @par Координатная система:
 *  - +X вправо, +Y вверх, +Z к наблюдателю (правая СК OpenGL).
 *  - Начало координат находится в центре масс системы тел
 *    после вызова BodySystem::transPointToSystem().
 *
 * @par Управление (клавиши, назначаются в Control):
 *  | Клавиша      | Действие                           |
 *  |--------------|------------------------------------|
 *  | W / S        | Движение вперёд / назад             |
 *  | A / D        | Стрейф влево / вправо               |
 *  | E / Q        | Движение вверх / вниз               |
 *  | Мышь         | Вращение (yaw, pitch)               |
 *  | Колёсико     | Изменение FOV (зум)                 |
 *  | Shift + ± | Изменение cameraMoveScale           |
 */
typedef struct Camera {

    glm::vec3 pos   = glm::vec3(0.0f, 0.0f, 120.0f); ///< Позиция камеры в мировых координатах [AU].
    glm::vec3 front = glm::vec3(0.0f, 0.0f, -1.0f);  ///< Единичный вектор направления взгляда.
    glm::vec3 up    = glm::vec3(0.0f, 1.0f,  0.0f);  ///< Единичный вектор «вверх» камеры.

    float yaw   = -90.0f; ///< Горизонтальный угол поворота [градусы]. -90° = смотрим в сторону -Z.
    float pitch =   0.0f; ///< Вертикальный угол поворота [градусы]. Ограничен ±89° во избежание gimbal lock.

    float lastX = 400.0f; ///< Предыдущая X-координата курсора мыши [пикс], для вычисления delta.
    float lastY = 300.0f; ///< Предыдущая Y-координата курсора мыши [пикс], для вычисления delta.

    float fov = 45.0f; ///< Вертикальный угол обзора [градусы]. Изменяется колёсиком.

    /**
     * @brief Вычисляет матрицу вида (View matrix).
     *
     * Строится как glm::lookAt(pos, pos + front, up).
     * Используется в Renderer для вычисления матрицы MVP каждого тела.
     *
     * @return Матрица вида 4×4 (float).
     */
    glm::mat4 getViewMatrix() const {
        return glm::lookAt(pos, pos + front, up);
    }

} Camera;

#endif // CAMERAINC
