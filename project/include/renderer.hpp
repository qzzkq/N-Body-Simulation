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
    
    GLFWwindow* window_ = nullptr; 

    // Инициализаторы геометрии
    void initSphereGeometry();
    void initCubeGeometry();  
    void initPointGeometry();
    void initTrailVAO_VBO(); 

    bool initWindow(int width, 
        int height, 
        const char* title, 
        bool fullscreen, 
        bool maximized);
    void initProgram(); 

    GLuint program_;
    GLint uModel_, uView_, uProj_, uColor_, uUseInstancing_;
    
    GLuint trailProgram_;
    GLuint pointProgram_;

    RenderMode mode_ = RenderMode::Cubes; 

    // VAO и VBO для сферы 
    GLuint sphereVAO_ = 0;
    GLuint sphereVBO_ = 0;
    size_t sphereVertexCount_ = 0;

    // VAO и VBO для куба 
    GLuint cubeVAO_ = 0;
    GLuint cubeVBO_ = 0;

    // VAO и VBO для точки 
    GLuint pointVAO_ = 0;
    GLuint pointVBO_ = 0;
    GLuint pointColorVBO_ = 0;

    // VBO для инстансинга (матрица модели + цвет)
    GLuint instanceModelVBO_ = 0;
    GLuint instanceColorVBO_ = 0;

    //VAO и VBO для трэйлов
    GLuint trailVAO_ = 0;
    GLuint trailVBO_ = 0;
    GLuint trailColorVBO_ = 0;

    // Выделенные размеры GPU-буферов — у каждого VBO своя ёмкость
    GLsizeiptr instanceModelCap_ = 0;
    GLsizeiptr instanceColorCap_ = 0;
    GLsizeiptr pointPosCap_      = 0;
    GLsizeiptr pointColorCap_    = 0;
    GLsizeiptr trailPosCap_      = 0;
    GLsizeiptr trailColorCap_    = 0;

    glm::mat4 viewMatrix_;
    glm::mat4 projectionMatrix_;
    float farPlane_ = 1.0e10f;

    bool successInit_ = false;
};
