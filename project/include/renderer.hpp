/**
 * @file renderer.hpp
 * @brief OpenGL-рендерер для 3D-визуализации N-body симуляции.
 *
 * Renderer инкапсулирует всё взаимодействие с OpenGL и GLFW:
 *  - Создание окна и контекста OpenGL 4.x.
 *  - Компиляцию GLSL-шейдеров (хранятся inline в renderer.cpp).
 *  - Рендеринг тел в трёх режимах: сферы, кубы, точки.
 *  - Рендеринг хвостов траекторий (GL_LINE_STRIP).
 *  - Инстансированный рендеринг (glDrawElementsInstanced) для N > 100.
 *
 * Логарифмический буфер глубины (реализован в шейдерах) позволяет
 * одновременно видеть объекты на расстоянии от ~1 AU до ~10¹⁰ AU.
 */

#pragma once
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <GLFW/glfw3.h>
#include <vector>
#include "object.hpp"
#include "graphic_state.hpp"
#include "camera.hpp"

/**
 * @enum RenderMode
 * @brief Режим геометрии для отрисовки тел.
 */
enum class RenderMode {
    Sphere, ///< Тела рисуются как сферы (UV-сфера, задаётся SPHERE_SUBDIVISIONS в renderer.cpp).
    Points, ///< Тела рисуются как GL_POINTS — максимальная производительность при N >> 10 000.
    Cubes   ///< Тела рисуются как кубы (отладочный/артистический режим).
};

/**
 * @class Renderer
 * @brief Управляет окном GLFW, OpenGL-контекстом и отрисовкой всей сцены.
 *
 * Типичный порядок использования:
 * @code
 * Renderer renderer;
 * renderer.init(1280, 720, "N-Body", false, true);
 * renderer.setProjection(65.0f, 1.0f, 1.0e10f);
 * renderer.setRenderMode(RenderMode::Sphere);
 *
 * while (!glfwWindowShouldClose(renderer.getWindow())) {
 *     // ... физика ...
 *     renderer.renderFrame(objs, graphics, cam);
 *     glfwPollEvents();
 * }
 * @endcode
 *
 * @note Renderer не является копируемым. Деструктор освобождает все GPU-ресурсы
 *       и вызывает glfwTerminate().
 */
class Renderer {
public:
    Renderer();
    ~Renderer();

    /**
     * @brief Инициализирует окно GLFW и OpenGL-контекст.
     *
     * Создаёт окно, инициализирует GLEW, компилирует шейдерные программы,
     * создаёт VAO/VBO для всех режимов геометрии и трейлов.
     * Устанавливает GLFW-коллбэк изменения размера окна.
     *
     * @param width       Ширина окна [пикс].
     * @param height      Высота окна [пикс].
     * @param title       Заголовок окна (отображается в ОС).
     * @param fullscreen  @c true — полноэкранный режим.
     * @param maximized   @c true — развернуть во весь экран при старте.
     * @return @c true при успехе, @c false при ошибке GLFW/GLEW.
     */
    bool init(int width, int height, const char* title, bool fullscreen, bool maximized);

    /**
     * @brief Устанавливает матрицу проекции (перспективная).
     *
     * Использует glm::perspective(). Матрица сохраняется в @c projectionMatrix_
     * и передаётся в шейдер при каждом вызове renderFrame().
     *
     * @param fov_deg  Вертикальный угол обзора [градусы], типично 45–65.
     * @param znear    Ближняя плоскость отсечения [AU], типично 1.0.
     * @param zfar     Дальняя плоскость отсечения [AU], типично 1e10.
     */
    void setProjection(float fov_deg, float znear, float zfar);

    /**
     * @brief Обновляет матрицу вида из текущего состояния камеры.
     *
     * @deprecated Используйте renderFrame() — он вызывает это внутри.
     * @param camera Состояние камеры.
     */
    void updateView(const Camera& camera);

    /**
     * @brief Выполняет полный цикл отрисовки одного кадра.
     *
     * Порядок вызовов:
     *  1. glClear(COLOR | DEPTH)
     *  2. updateView(cam)
     *  3. drawObjects(objs, graphics, cam)
     *  4. drawTrails(objs, graphics, cam)
     *  5. glfwSwapBuffers()
     *
     * @param objs     Массив физических тел (позиции, радиусы).
     * @param graphics Массив визуальных состояний (цвета, хвосты).
     * @param cam      Текущее состояние камеры.
     */
    void renderFrame(const std::vector<Object>& objs,
                     const std::vector<GraphicState>& graphics,
                     const Camera& cam);

    /**
     * @brief Отрисовывает все тела в текущем RenderMode через инстансинг.
     *
     * Заполняет instanceModelVBO_ (матрицы модели 4×4) и instanceColorVBO_
     * из обновлённых позиций тел, затем вызывает glDrawElementsInstanced
     * (Sphere/Cubes) или glDrawArraysInstanced (Points).
     *
     * @param objs     Массив тел.
     * @param graphics Массив визуальных состояний.
     * @param cam      Камера (для расчёта расстояния LOD при Points).
     */
    void drawObjects(const std::vector<Object>& objs,
                     const std::vector<GraphicState>& graphics,
                     const Camera& cam);

    /**
     * @brief Отрисовывает хвосты траекторий всех тел.
     *
     * Для каждого тела с непустым GraphicState::trail рисует GL_LINE_STRIP.
     * Цвет совпадает с цветом тела, альфа плавно убывает к «хвосту».
     *
     * @param objs     Массив тел.
     * @param graphics Массив визуальных состояний (используется trail).
     * @param cam      Камера.
     */
    void drawTrails(const std::vector<Object>& objs,
                    const std::vector<GraphicState>& graphics,
                    const Camera& cam);

    /**
     * @brief Устанавливает режим геометрии тел.
     * @param mode Новый режим (Sphere / Points / Cubes).
     */
    void setRenderMode(RenderMode mode) { mode_ = mode; }

    /**
     * @brief Возвращает указатель на окно GLFW.
     *
     * Используется в main.cpp и Control для регистрации коллбэков,
     * опроса состояния клавиш и swap buffers.
     * @return Указатель на GLFWwindow. Не владеет.
     */
    GLFWwindow* getWindow() const { return window_; }

    /**
     * @brief Обрабатывает изменение размера окна.
     *
     * Вызывается из GLFW framebuffer size callback.
     * Обновляет glViewport и пересчитывает матрицу проекции
     * с сохранённым FOV и плоскостями отсечения.
     *
     * @param w Новая ширина [пикс].
     * @param h Новая высота [пикс].
     */
    void resizeWindow(int w, int h);

private:
    /**
     * @brief Компилирует и линкует GLSL-программу из двух шейдеров.
     * @param vs Исходный код вершинного шейдера.
     * @param fs Исходный код фрагментного шейдера.
     * @return OpenGL-хэндл программы (GLuint).
     */
    GLuint compileProgram(const char* vs, const char* fs);

    /**
     * @brief Настраивает attribute pointers для инстансированного рендеринга.
     *
     * Подключает instanceModelVBO_ (mat4, layout 4–7, divisor=1) и
     * instanceColorVBO_ (vec4, layout 8, divisor=1) к заданному VAO.
     * @param vao Целевой VAO.
     */
    void setupInstanceAttribs(GLuint vao);

    GLFWwindow* window_ = nullptr; ///< Окно GLFW и OpenGL-контекст.

    /// @name Инициализаторы геометрии
    /// @{
    void initSphereGeometry(); ///< Генерирует UV-сферу, загружает в sphereVAO_/sphereVBO_.
    void initCubeGeometry();   ///< Генерирует единичный куб, загружает в cubeVAO_/cubeVBO_.
    void initPointGeometry();  ///< Создаёт VAO для режима точек.
    void initTrailVAO_VBO();   ///< Создаёт VAO/VBO для рисования хвостов (динамический буфер).
    /// @}

    bool initWindow(int width, int height, const char* title,
                    bool fullscreen, bool maximized); ///< Вспомогательный метод создания окна GLFW.
    void initProgram(); ///< Компилирует все шейдерные программы (основная, трейлы, точки).

    // ── Шейдерные программы ──────────────────────────────────────────────
    GLuint program_;       ///< Основная программа: Sphere и Cubes (инстансинг).
    GLuint trailProgram_;  ///< Программа для хвостов траекторий (GL_LINE_STRIP).
    GLuint pointProgram_;  ///< Программа для режима Points (GL_POINTS с gl_PointSize).

    // ── Uniforms основной программы ──────────────────────────────────────
    GLint uModel_;          ///< Локация uniform mat4 uModel.
    GLint uView_;           ///< Локация uniform mat4 uView.
    GLint uProj_;           ///< Локация uniform mat4 uProj.
    GLint uColor_;          ///< Локация uniform vec4 uColor.
    GLint uUseInstancing_;  ///< Локация uniform bool uUseInstancing.

    RenderMode mode_ = RenderMode::Cubes; ///< Текущий режим отрисовки тел.

    // ── VAO / VBO геометрий ──────────────────────────────────────────────
    GLuint sphereVAO_ = 0;       ///< VAO сферы.
    GLuint sphereVBO_ = 0;       ///< VBO сферы (вершины + нормали).
    size_t sphereVertexCount_ = 0; ///< Количество вершин сферы.

    GLuint cubeVAO_ = 0;         ///< VAO куба.
    GLuint cubeVBO_ = 0;         ///< VBO куба.

    GLuint pointVAO_ = 0;        ///< VAO для режима точек.
    GLuint pointVBO_ = 0;        ///< VBO позиций точек (обновляется каждый кадр).
    GLuint pointColorVBO_ = 0;   ///< VBO цветов точек (обновляется каждый кадр).

    // ── VBO инстансинга ──────────────────────────────────────────────────
    GLuint instanceModelVBO_ = 0; ///< VBO матриц модели (4×4 float, по одной на тело).
    GLuint instanceColorVBO_ = 0; ///< VBO цветов тел (vec4 float, по одному на тело).

    // ── Хвосты ───────────────────────────────────────────────────────────
    GLuint trailVAO_ = 0; ///< VAO для трейлов (переиспользуется для всех тел).
    GLuint trailVBO_ = 0; ///< VBO позиций точек хвоста (динамический, GL_DYNAMIC_DRAW).

    // ── Матрицы ──────────────────────────────────────────────────────────
    glm::mat4 viewMatrix_;       ///< Текущая матрица вида (обновляется каждый кадр).
    glm::mat4 projectionMatrix_; ///< Матрица проекции (обновляется при ресайзе).
    float farPlane_ = 1.0e10f;   ///< Дальняя плоскость (хранится для пересчёта при ресайзе).

    bool successInit_ = false;   ///< @c true если init() завершился успешно.
};
