#include "renderer.hpp"
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>
#include <cmath>

namespace {
    constexpr double kG = 6.6743e-11;
}

Renderer::Renderer(int /*w*/, int /*h*/) { // конструктор renderer'а 
    this->program_ = compileProgram();
    glUseProgram(program_);
    this->uModel_ = glGetUniformLocation(program_, "model");
    this->uView_  = glGetUniformLocation(program_, "view");
    this->uProj_  = glGetUniformLocation(program_, "projection");
    this->uColor_ = glGetUniformLocation(program_, "objectColor");

    initCube(this->cubeMesh_); 
    initPoint(this->pointMesh_);
}

Renderer::~Renderer() { // деструктор - удаление шейдерной программы
    if (this->program_) 
        glDeleteProgram(program_);

    if (this->cubeMesh_.VAO) {
        glDeleteVertexArrays(1, &cubeMesh_.VAO);
        glDeleteBuffers(1, &cubeMesh_.VBO);
    }
    if (this->pointMesh_.VAO) {
        glDeleteVertexArrays(1, &pointMesh_.VAO);
        glDeleteBuffers(1, &pointMesh_.VBO);
    }
}

void Renderer::setProjection(float fov_deg, // метод, устанавливающий параметры обзора камеры 
                             float aspect,
                             float znear,
                             float zfar)
{
    glUseProgram(program_);
    glm::mat4 P = glm::perspective(glm::radians(fov_deg), aspect, znear, zfar);
    glUniformMatrix4fv(uProj_, 1, GL_FALSE, glm::value_ptr(P));
}

void Renderer::updateView(const glm::vec3& pos, // изменение обзора 
                          const glm::vec3& front,
                          const glm::vec3& up)
{
    glUseProgram(program_);
    glm::mat4 V = glm::lookAt(pos, pos + front, up);
    glUniformMatrix4fv(uView_, 1, GL_FALSE, glm::value_ptr(V));
}

void Renderer::drawObjects(const std::vector<Object>& objs) const { // функция отрисовки объекта 
    glUseProgram(program_);

    const Mesh* currentMesh = nullptr;
    GLenum drawMode = GL_TRIANGLES;

    switch (this->mode_) {
        case Mode::Points: 
            currentMesh = &pointMesh_;
            drawMode = GL_POINTS;
            glEnable(GL_PROGRAM_POINT_SIZE); // Включаем возможность менять размер точки
            glPointSize(6.0f); // Устанавливаем размер точки
            break;
        case Mode::Cubes: 
        default: 
            currentMesh = &cubeMesh_;
            drawMode = GL_TRIANGLES;
            break;
    }

    if (currentMesh && currentMesh->VAO) {
        glBindVertexArray(currentMesh->VAO);
        
        for (const auto& obj : objs) {
            float scale = static_cast<float>(obj.GetRadius());
            float halfSide = glm::max(scale, 0.0001f); 
            glm::mat4 M(1.0f);
            M = glm::translate(M, glm::vec3(obj.GetPos()));
            if (mode_ == Mode::Cubes) {
                M = glm::scale(M, glm::vec3(halfSide));
            }
            glUniformMatrix4fv(uModel_, 1, GL_FALSE, glm::value_ptr(M));
            glUniform4f(uColor_, obj.GetColor().r, obj.GetColor().g, obj.GetColor().b, obj.GetColor().a);
            glDrawArrays(drawMode, 0, currentMesh->vertexCount);
        }
        glBindVertexArray(0); // Отвязываем VAO 
        if (mode_ == Mode::Points) {
            glDisable(GL_PROGRAM_POINT_SIZE);
        }
    }
}

GLuint Renderer::compileProgram() { // компиляция шейдерной программы для GPU
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
    GLuint v = compile(GL_VERTEX_SHADER,   this->vs); // компиляция вершинного шейдера
    GLuint f = compile(GL_FRAGMENT_SHADER, this->fs); // компиляция фрагментного шейдера
    GLuint p = glCreateProgram();
    glAttachShader(p, v);
    glAttachShader(p, f);
    glLinkProgram(p);
    glDeleteShader(v);
    glDeleteShader(f);
    return p;
}

void Renderer::initCube(Mesh& mesh) {
    // Вектор вершин для куба 2x2x2 (от -1 до 1) - 36 вершин (12 треугольников)
    std::vector<float> vertices = {
        //   X      Y      Z
        // back
        -1.f,-1.f,-1.f, 1.f, 1.f,-1.f, 1.f,-1.f,-1.f,
        -1.f,-1.f,-1.f, -1.f, 1.f,-1.f, 1.f, 1.f,-1.f,
        // front
        -1.f,-1.f, 1.f, 1.f,-1.f, 1.f, 1.f, 1.f, 1.f,
        -1.f,-1.f, 1.f, 1.f, 1.f, 1.f, -1.f, 1.f, 1.f,
        // left
        -1.f,-1.f,-1.f, -1.f,-1.f, 1.f, -1.f, 1.f, 1.f,
        -1.f,-1.f,-1.f, -1.f, 1.f, 1.f, -1.f, 1.f,-1.f,
        // right
         1.f,-1.f,-1.f, 1.f, 1.f, 1.f, 1.f,-1.f, 1.f,
         1.f,-1.f,-1.f, 1.f, 1.f,-1.f, 1.f, 1.f, 1.f,
        // bottom
        -1.f,-1.f,-1.f, 1.f,-1.f, 1.f, 1.f,-1.f,-1.f,
        -1.f,-1.f,-1.f, -1.f,-1.f, 1.f, 1.f,-1.f, 1.f,
        // top
        -1.f, 1.f,-1.f, 1.f, 1.f,-1.f, 1.f, 1.f, 1.f,
        -1.f, 1.f,-1.f, 1.f, 1.f, 1.f, -1.f, 1.f, 1.f
    };
    initMesh(mesh, vertices, 3);
}

void Renderer::initPoint(Mesh& mesh) {
    std::vector<float> vertices = { 0.0f, 0.0f, 0.0f };
    initMesh(mesh, vertices, 3);
}

void Renderer::initMesh(Mesh& mesh, const std::vector<float>& vertices, int vertexComponents) { // функция, которая загружает данные для визуализации
    if (vertices.empty()) 
        return;
    if (mesh.VAO) {
        glDeleteVertexArrays(1, &(mesh.VAO));
        glDeleteBuffers(1, &(mesh.VBO));
    }
    mesh.vertexCount = static_cast<GLsizei>(vertices.size() / vertexComponents); 
    glGenVertexArrays(1, &mesh.VAO);
    glGenBuffers(1, &mesh.VBO);
    glBindVertexArray(mesh.VAO);
    glBindBuffer(GL_ARRAY_BUFFER, mesh.VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, vertexComponents, GL_FLOAT, GL_FALSE, vertexComponents * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0); 
}