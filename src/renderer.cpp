#include "renderer.hpp"
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>
#include <cmath>

namespace {
    constexpr double kG = 6.6743e-11;
}

Renderer::Renderer(int /*w*/, int /*h*/, const char* vs, const char* fs) {
    program_ = compileProgram(vs, fs);
    glUseProgram(program_);
    uModel_ = glGetUniformLocation(program_, "model");
    uView_  = glGetUniformLocation(program_, "view");
    uProj_  = glGetUniformLocation(program_, "projection");
    uColor_ = glGetUniformLocation(program_, "objectColor");

    glGenVertexArrays(1, &gridVAO_);
    glGenBuffers(1, &gridVBO_);
    glBindVertexArray(gridVAO_);
    glBindBuffer(GL_ARRAY_BUFFER, gridVBO_);
    glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);
}

Renderer::~Renderer() {
    if (gridVBO_) glDeleteBuffers(1, &gridVBO_);
    if (gridVAO_) glDeleteVertexArrays(1, &gridVAO_);
    if (program_) glDeleteProgram(program_);
}

void Renderer::setProjection(float fov_deg,
                             float aspect,
                             float znear,
                             float zfar)
{
    glUseProgram(program_);
    glm::mat4 P = glm::perspective(glm::radians(fov_deg), aspect, znear, zfar);
    glUniformMatrix4fv(uProj_, 1, GL_FALSE, glm::value_ptr(P));
}

void Renderer::updateView(const glm::vec3& pos,
                          const glm::vec3& front,
                          const glm::vec3& up)
{
    glUseProgram(program_);
    glm::mat4 V = glm::lookAt(pos, pos + front, up);
    glUniformMatrix4fv(uView_, 1, GL_FALSE, glm::value_ptr(V));
}

void Renderer::updateGrid(float size, int divisions, const std::vector<Object>& objs) {
    auto verts = createGridVertices(size, divisions, objs);
    gridVertexCount_ = verts.size() / 3;

    glBindBuffer(GL_ARRAY_BUFFER, gridVBO_);
    glBufferData(GL_ARRAY_BUFFER, verts.size()*sizeof(float), verts.data(), GL_DYNAMIC_DRAW);
}

void Renderer::drawGrid() const {
    glUseProgram(program_);
    glm::mat4 M(1.0f);
    glUniformMatrix4fv(uModel_, 1, GL_FALSE, glm::value_ptr(M));
    glUniform4f(uColor_, 1.0f, 1.0f, 1.0f, 0.25f);

    glBindVertexArray(gridVAO_);
    glDrawArrays(GL_LINES, 0, static_cast<GLint>(gridVertexCount_));
    glBindVertexArray(0);
}

// ===== вот тут объединены 3 режима =====
void Renderer::drawObjects(const std::vector<Object>& objs) const {
    glUseProgram(program_);

    switch (mode_) {
    case Mode::Points: {
        glPointSize(6.0f);
        for (const auto& obj : objs) {
            glm::mat4 M(1.0f);
            M = glm::translate(M, glm::vec3(obj.position));
            glUniformMatrix4fv(uModel_, 1, GL_FALSE, glm::value_ptr(M));
            glUniform4f(uColor_, obj.color.r, obj.color.g, obj.color.b, obj.color.a);

            glBindVertexArray(obj.VAO);
            glDrawArrays(GL_POINTS, 0, 1);
        }
        break;
    }

    case Mode::Cubes: {
        // рисуем один юнит-куб и масштабируем под радиус
        static GLuint cubeVAO = 0, cubeVBO = 0;
        if (cubeVAO == 0) {
            float cubeVerts[] = {
                // back
                -1.f,-1.f,-1.f,  1.f, 1.f,-1.f,  1.f,-1.f,-1.f,
                -1.f,-1.f,-1.f, -1.f, 1.f,-1.f,  1.f, 1.f,-1.f,
                // front
                -1.f,-1.f, 1.f,  1.f,-1.f, 1.f,  1.f, 1.f, 1.f,
                -1.f,-1.f, 1.f,  1.f, 1.f, 1.f, -1.f, 1.f, 1.f,
                // left
                -1.f,-1.f,-1.f, -1.f,-1.f, 1.f, -1.f, 1.f, 1.f,
                -1.f,-1.f,-1.f, -1.f, 1.f, 1.f, -1.f, 1.f,-1.f,
                // right
                 1.f,-1.f,-1.f,  1.f, 1.f, 1.f,  1.f,-1.f, 1.f,
                 1.f,-1.f,-1.f,  1.f, 1.f,-1.f,  1.f, 1.f, 1.f,
                // bottom
                -1.f,-1.f,-1.f,  1.f,-1.f, 1.f,  1.f,-1.f,-1.f,
                -1.f,-1.f,-1.f, -1.f,-1.f, 1.f,  1.f,-1.f, 1.f,
                // top
                -1.f, 1.f,-1.f,  1.f, 1.f,-1.f,  1.f, 1.f, 1.f,
                -1.f, 1.f,-1.f,  1.f, 1.f, 1.f, -1.f, 1.f, 1.f
            };
            glGenVertexArrays(1, &cubeVAO);
            glGenBuffers(1, &cubeVBO);
            glBindVertexArray(cubeVAO);
            glBindBuffer(GL_ARRAY_BUFFER, cubeVBO);
            glBufferData(GL_ARRAY_BUFFER, sizeof(cubeVerts), cubeVerts, GL_STATIC_DRAW);
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(float), (void*)0);
            glEnableVertexAttribArray(0);
            glBindVertexArray(0);
        }

        for (const auto& obj : objs) {
            float halfSide = glm::max(static_cast<float>(obj.radius), 0.0001f);

            glm::mat4 M(1.0f);
            M = glm::translate(M, glm::vec3(obj.position));
            M = glm::scale(M, glm::vec3(halfSide));
            glUniformMatrix4fv(uModel_, 1, GL_FALSE, glm::value_ptr(M));
            glUniform4f(uColor_, obj.color.r, obj.color.g, obj.color.b, obj.color.a);

            glBindVertexArray(cubeVAO);
            glDrawArrays(GL_TRIANGLES, 0, 36);
        }
        break;
    }

    case Mode::Sphere:
    default: {
        // твой исходный вариант — просто рисуем VAO объекта
        for (const auto& obj : objs) {
            glm::mat4 M(1.0f);
            M = glm::translate(M, glm::vec3(obj.position));
            glUniformMatrix4fv(uModel_, 1, GL_FALSE, glm::value_ptr(M));
            glUniform4f(uColor_, obj.color.r, obj.color.g, obj.color.b, obj.color.a);

            glBindVertexArray(obj.VAO);
            glDrawArrays(GL_TRIANGLES, 0, static_cast<GLint>(obj.vertexCount / 3));
        }
        break;
    }
    }
}

GLuint Renderer::compileProgram(const char* vs, const char* fs) {
    auto compile = [](GLenum type, const char* src){
        GLuint s = glCreateShader(type);
        glShaderSource(s, 1, &src, nullptr);
        glCompileShader(s);
        GLint ok=0; glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
        if(!ok){
            char log[1024];
            glGetShaderInfoLog(s,1024,nullptr,log);
            std::cerr << "Shader compile error: " << log << std::endl;
        }
        return s;
    };
    GLuint v = compile(GL_VERTEX_SHADER,   vs);
    GLuint f = compile(GL_FRAGMENT_SHADER, fs);
    GLuint p = glCreateProgram();
    glAttachShader(p, v);
    glAttachShader(p, f);
    glLinkProgram(p);
    glDeleteShader(v);
    glDeleteShader(f);
    return p;
}

std::vector<float> Renderer::createGridVertices(float size, int divisions,
                                                const std::vector<Object>& objs) {
    std::vector<float> vertices;
    float step = size / divisions;
    float half = size * 0.5f;

    // ... твой код генерации сетки и смещения (оставь как был) ...
    // (я его не переписываю тут целиком, у тебя он уже есть в файле)
    // главное — эта функция должна остаться такой же, как была у тебя.
    // ↓↓↓
    // (скопируй сюда свой имеющийся createGridVertices из проекта)
    return vertices;
}
