#pragma once
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp> 
#include <vector>
#include "object.hpp"

typedef struct Mesh {
    GLuint VAO = 0;
    GLuint VBO = 0;
    GLsizei vertexCount = 0;
} Mesh; 

class Renderer {
public:
    // режим визуализации, выбираем его из main
    enum class Mode {
        Points,      // одна точка на объект
        Cubes        // куб вместо объекта
    };

    Renderer(int width, int height);
    ~Renderer();

    void setProjection(float fov_deg, float aspect, float znear, float zfar);
    void updateView(const glm::vec3& pos,
                    const glm::vec3& front,
                    const glm::vec3& up);

    void drawObjects(const std::vector<Object>& objs) const;
    void setMode(Mode m) { mode_ = m; }
    Mode getMode() const { return mode_; }

private:
    GLuint program_ = 0;
    GLint  uModel_ = -1, uView_ = -1, uProj_ = -1, uColor_ = -1;

    // Вершинный шейдер
    const char* vs = R"glsl(
        #version 330 core
        layout(location=0) in vec3 aPos;
        uniform mat4 model;
        uniform mat4 view;
        uniform mat4 projection;
        void main() {
            gl_Position = projection * view * model * vec4(aPos, 1.0);
        }
    )glsl";
    
    // Фрагментный шейдер
    const char* fs = R"glsl(
        #version 330 core
        out vec4 FragColor;
        uniform vec4 objectColor;
        void main() {
            FragColor = objectColor;
        }
    )glsl";

    struct Mesh pointMesh_; 
    struct Mesh sphereMesh_; 
    struct Mesh cubeMesh_; 

    void initMesh(Mesh& mesh, const std::vector<float>& vertices, int vertexComponents); // метод загрузки данных в OpenGL
    void initCube(Mesh& cubeMesh_); // метод инициализации кастомного примитива "Куб"
    void initPoint(Mesh& pointMesh_); // метод инициализации кастомного примитика "точка"

    // выбранный режим отрисовки
    Mode mode_ = Mode::Points;

    GLuint compileProgram();
};
