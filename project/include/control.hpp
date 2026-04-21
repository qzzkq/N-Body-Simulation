/**
 * @file control.hpp
 * @brief Обработчик пользовательского ввода: клавиатура, мышь, создание тел.
 *
 * Control регистрирует GLFW-коллбэки и транслирует события ввода
 * в изменения SimState (пауза, масштаб времени) и Camera (движение, вращение).
 *
 * Использует паттерн «статическая обёртка + user pointer»:
 * GLFW требует статические функции для коллбэков, поэтому
 * каждый коллбэк достаёт указатель на экземпляр Control через
 * glfwGetWindowUserPointer() и делегирует вызов instance-методу.
 */

#pragma once
struct GLFWwindow;
#include <glm/glm.hpp>
#include <vector>
#include "object.hpp"
#include "camera.hpp"
#include "state.hpp"

/**
 * @class Control
 * @brief Управляет вводом пользователя и синхронизирует состояние камеры и симуляции.
 *
 * Создаётся в main() после инициализации Renderer.
 * Хранит ссылки — не владеет объектами.
 *
 * @par Полная таблица управления:
 *  | Клавиша / событие        | Действие                                  |
 *  |--------------------------|-------------------------------------------|
 *  | W / A / S / D            | Движение камеры (вперёд/влево/назад/вправо)|
 *  | E / Q                    | Движение камеры вверх / вниз              |
 *  | Мышь (движение)          | Вращение камеры (yaw / pitch)             |
 *  | Колёсико                 | Изменение FOV                             |
 *  | Space                    | Пауза / возобновление симуляции           |
 *  | Escape                   | Завершение работы (state.running = false) |
 *  | + / - (без Shift)        | Изменение timeScale (×1.5 / ÷1.5)        |
 *  | Shift + = / Shift + -    | Увеличить / уменьшить cameraMoveScale     |
 *  | Shift (удержание)        | «Заморозка» времени без паузы             |
 *
 * @note Режим мыши: GLFW_CURSOR_DISABLED (скрыт, захвачен окном).
 *       Это означает что для выхода из приложения нужно нажать Escape.
 */
class Control {
public:
    /**
     * @brief Конструктор. Принимает ссылки на все необходимые объекты состояния.
     *
     * @param window  Указатель на окно GLFW. Используется для регистрации коллбэков.
     * @param objs    Массив тел симуляции (для будущего интерактивного добавления тел).
     * @param cam     Камера наблюдателя (изменяется при вводе).
     * @param state   Состояние симуляции (pause, timeScale, running).
     */
    Control(GLFWwindow* window,
            std::vector<Object>& objs,
            Camera& cam,
            SimState& state);

    /**
     * @brief Регистрирует GLFW-коллбэки и захватывает курсор.
     *
     * Вызывает:
     *  - glfwSetWindowUserPointer(window_, this)
     *  - glfwSetKeyCallback / ScrollCallback / CursorPosCallback
     *  - glfwSetInputMode(GLFW_CURSOR, GLFW_CURSOR_DISABLED)
     *
     * Должен быть вызван один раз после создания экземпляра.
     */
    void attach();

    /**
     * @brief Обновляет позицию камеры по нажатым клавишам.
     *
     * Опрашивает состояние клавиш W/A/S/D/E/Q через glfwGetKey()
     * и перемещает Camera::pos вдоль front/right/up векторов.
     * Скорость движения = cameraMoveScale_ * state_.deltaTime.
     *
     * Вызывается **каждый кадр** из главного цикла, в отличие от
     * onKey() (который реагирует на события нажатия/отпускания).
     */
    void updateCameraFromKeys();

    /**
     * @brief Возвращает текущий множитель скорости движения камеры.
     *
     * Изменяется клавишами Shift+= / Shift+-. Отображается в заголовке окна.
     * @return Множитель cameraMoveScale [безразмерный, > 0].
     */
    float getCameraMoveScale() const { return cameraMoveScale_; }

private:
    /**
     * @brief Статическая обёртка GLFW для обработки клавиш.
     * Делегирует вызов в Control::onKey() через glfwGetWindowUserPointer.
     */
    static void KeyCB(GLFWwindow* w, int key, int scancode, int action, int mods);

    /**
     * @brief Статическая обёртка GLFW для обработки колёсика мыши.
     * Делегирует вызов в Control::onScroll().
     */
    static void ScrollCB(GLFWwindow* w, double xoffset, double yoffset);

    /**
     * @brief Статическая обёртка GLFW для обработки движения мыши.
     * Делегирует вызов в Control::onCursorPos().
     */
    static void CursorPosCB(GLFWwindow* w, double xpos, double ypos);

    /**
     * @brief Обрабатывает события клавиш: пауза, выход, timeScale, cameraMoveScale.
     * @param key      GLFW-код клавиши (GLFW_KEY_*).
     * @param scancode Аппаратный скан-код (не используется).
     * @param action   GLFW_PRESS / GLFW_RELEASE / GLFW_REPEAT.
     * @param mods     Битовая маска модификаторов (GLFW_MOD_SHIFT и др.).
     */
    void onKey(int key, int scancode, int action, int mods);

    /**
     * @brief Обрабатывает прокрутку колёсика: изменяет Camera::fov.
     * @param xoffset Горизонтальная прокрутка (обычно 0).
     * @param yoffset Вертикальная прокрутка (±1 за щелчок).
     */
    void onScroll(double xoffset, double yoffset);

    /**
     * @brief Обрабатывает движение курсора мыши: вращает камеру (yaw, pitch).
     *
     * Вычисляет delta от Camera::lastX/lastY, умножает на sensitivity
     * и обновляет Camera::yaw / Camera::pitch. Затем пересчитывает
     * Camera::front через тригонометрию.
     *
     * @param xpos  Текущая X-координата курсора [пикс].
     * @param ypos  Текущая Y-координата курсора [пикс].
     */
    void onCursorPos(double xpos, double ypos);

    GLFWwindow*          window_;          ///< Окно GLFW (не владеет).
    std::vector<Object>& objs_;            ///< Массив тел (для интерактивного добавления).
    Camera&              camera_;          ///< Камера наблюдателя.
    SimState&            state_;           ///< Состояние симуляции.
    float                cameraMoveScale_; ///< Множитель скорости WASD. Регулируется Shift+±.
};
