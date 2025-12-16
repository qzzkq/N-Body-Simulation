#include "renderer.hpp"
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>
#include <cmath>

namespace {
    constexpr double kG = 6.6743e-11;
}

// перевод сферических координат в декартовы 
static glm::vec3 sphericalToCartesian(float r, float theta, float phi) {
    return glm::vec3(
        r * sin(theta) * cos(phi),
        r * cos(theta),
        r * sin(theta) * sin(phi)
    );
}

//==PUBLIC BLOCK==

Renderer::Renderer(int /*w*/, int /*h*/, const char* vs, const char* fs) {
    program_ = compileProgram(vs, fs);
    glUseProgram(program_);
    
    uModel_ = glGetUniformLocation(program_, "model");
    uView_  = glGetUniformLocation(program_, "view");
    uProj_  = glGetUniformLocation(program_, "projection");
    uColor_ = glGetUniformLocation(program_, "objectColor");

    initSphereGeometry();
    initCubeGeometry();
    initPointGeometry();
}

Renderer::~Renderer() {
    glDeleteVertexArrays(1, &sphereVAO_); glDeleteBuffers(1, &sphereVBO_);
    glDeleteVertexArrays(1, &cubeVAO_);   glDeleteBuffers(1, &cubeVBO_);
    glDeleteVertexArrays(1, &pointVAO_);  glDeleteBuffers(1, &pointVBO_);
    glDeleteProgram(program_);
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

void Renderer::drawObjects(const std::vector<Object>& objs) const {
    glUseProgram(program_);

    // Выбираем отрисовку в зависимости от режима 
    if (mode_ == RenderMode::Points) {
        glPointSize(4.0f);
        glBindVertexArray(pointVAO_);
        for (const auto& obj : objs) {
             glm::mat4 M(1.0f);
             M = glm::translate(M, glm::vec3(obj.position));
             glUniformMatrix4fv(uModel_, 1, GL_FALSE, glm::value_ptr(M));
             glUniform4f(uColor_, obj.color.r, obj.color.g, obj.color.b, obj.color.a);
             glDrawArrays(GL_POINTS, 0, 1);
        }
        glBindVertexArray(0);
    }
    else if (mode_ == RenderMode::Cubes) {
        glBindVertexArray(cubeVAO_);
        for (const auto& obj : objs) {
            glm::mat4 M(1.0f);
            M = glm::translate(M, glm::vec3(obj.position));
            float scale = std::max(obj.radius, 0.001f);
            M = glm::scale(M, glm::vec3(scale));
            
            glUniformMatrix4fv(uModel_, 1, GL_FALSE, glm::value_ptr(M));
            glUniform4f(uColor_, obj.color.r, obj.color.g, obj.color.b, obj.color.a);
            glDrawArrays(GL_TRIANGLES, 0, 36);
        }
        glBindVertexArray(0);
    }
    else { 
        glBindVertexArray(sphereVAO_);
        for (const auto& obj : objs) {
            glm::mat4 M(1.0f);
            M = glm::translate(M, glm::vec3(obj.position));
            float scale = std::max(obj.radius, 0.1f);
            M = glm::scale(M, glm::vec3(scale));
            
            glUniformMatrix4fv(uModel_, 1, GL_FALSE, glm::value_ptr(M));
            glUniform4f(uColor_, obj.color.r, obj.color.g, obj.color.b, obj.color.a);
            glDrawArrays(GL_TRIANGLES, 0, (GLsizei)sphereVertexCount_);
        }
        glBindVertexArray(0);
    }
}

//==PRIVATE BLOCK==

// Private-методы для компиляции шефдерной программы 

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

// Private-методы для инициализации отрисовочных параметров 

void Renderer::initPointGeometry() {
    // Точка в центре (0,0,0)
    float pointVert[] = { 0.0f, 0.0f, 0.0f };
    glGenVertexArrays(1, &pointVAO_);
    glGenBuffers(1, &pointVBO_);
    glBindVertexArray(pointVAO_);
    glBindBuffer(GL_ARRAY_BUFFER, pointVBO_);
    glBufferData(GL_ARRAY_BUFFER, sizeof(pointVert), pointVert, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);
}

void Renderer::initCubeGeometry() {
    // Стандартный куб -1..1
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
    glGenVertexArrays(1, &cubeVAO_);
    glGenBuffers(1, &cubeVBO_);
    glBindVertexArray(cubeVAO_);
    glBindBuffer(GL_ARRAY_BUFFER, cubeVBO_);
    glBufferData(GL_ARRAY_BUFFER, sizeof(cubeVerts), cubeVerts, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);
}

void Renderer::initSphereGeometry() {
    std::vector<float> vertices;
    const int stacks = 12;  // Детализация
    const int sectors = 12;
    const float radius = 1.0f; // Единичная сфера

    for (int i = 0; i < stacks; ++i) {
        float theta1 = (float(i) / stacks) * glm::pi<float>();
        float theta2 = (float(i + 1) / stacks) * glm::pi<float>();

        for (int j = 0; j < sectors; ++j) {
            float phi1 = (float(j) / sectors) * 2.0f * glm::pi<float>();
            float phi2 = (float(j + 1) / sectors) * 2.0f * glm::pi<float>();

            glm::vec3 v1 = sphericalToCartesian(radius, theta1, phi1);
            glm::vec3 v2 = sphericalToCartesian(radius, theta1, phi2);
            glm::vec3 v3 = sphericalToCartesian(radius, theta2, phi1);
            glm::vec3 v4 = sphericalToCartesian(radius, theta2, phi2);

            // Два треугольника на сектор
            auto addV = [&](const glm::vec3& v) {
                vertices.push_back(v.x); vertices.push_back(v.y); vertices.push_back(v.z);
            };

            addV(v1); addV(v2); addV(v3);
            addV(v2); addV(v4); addV(v3);
        }
    }
    
    sphereVertexCount_ = vertices.size() / 3;

    glGenVertexArrays(1, &sphereVAO_);
    glGenBuffers(1, &sphereVBO_);
    
    glBindVertexArray(sphereVAO_);
    glBindBuffer(GL_ARRAY_BUFFER, sphereVBO_);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);
    
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);
}