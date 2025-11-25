#include "renderer.hpp"
#include "constants.hpp"
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>
#include <cmath>

const char* vertexShaderSource = R"glsl(
#version 330 core
layout(location=0) in vec3 aPos;
uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
void main() {
    gl_Position = projection * view * model * vec4(aPos, 1.0);
}
)glsl";

const char* fragmentShaderSource = R"glsl(
#version 330 core
out vec4 FragColor;
uniform vec4 objectColor;
void main() {
    FragColor = objectColor;
}
)glsl";

Renderer::Renderer(int w, int h) {
    program_ = compileProgram(vertexShaderSource, fragmentShaderSource);
    initShader();
    initPrimitives();
}

Renderer::~Renderer() { // cleaning GPU memory
    if (sphereVBO_) glDeleteBuffers(1, &sphereVBO_); 
    if (sphereVAO_) glDeleteVertexArrays(1, &sphereVAO_);
    if (cubeVBO_) glDeleteBuffers(1, &cubeVBO_); 
    if (cubeVAO_) glDeleteVertexArrays(1, &cubeVAO_);
    if (program_) glDeleteProgram(program_);
}

Renderer::compileProgram(const char* vs, const char* fs) {
    auto compile = [](GLenum type, const char* src){
        GLuint s = glCreateShader(type);
        glShaderSource(s, 1, &src, nullptr);
        glCompileShader(s);
        GLint ok = 0; 
        glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
        if(!ok){
            char log[1024];
            glGetShaderInfoLog(s, 1024, nullptr, log);
            std::cerr << "Shader compile error: " << log << std::endl;
        }
        return s;
    };
    GLuint verSh = compile(GL_VERTEX_SHADER, vs);
    GLuint fragSh = compile(GL_FRAGMENT_SHADER, fs);
    GLuint prog = glCreateProgram();
    glAttachShader(prog, verSh);
    glAttachShader(prog, fragSh);
    glLinkProgram(prog);
    glDeleteShader(verSh);
    glDeleteShader(fragSh);

    GLint success;
    glGetProgramiv(prog, GL_LINK_STATUS, &success);
    if (!success) {
        char infoLog[1024];
        glGetProgramInfoLog(prog, 1024, NULL, infoLog);
        std::cerr << "Shader program linking failed: " << infoLog << std::endl;
        return 0;
    }
    return prog;
}

void Renderer::initShader() {
    glUseProgram(program_);
    uModel_ = glGetUniformLocation(program_, "model");
    uView_  = glGetUniformLocation(program_, "view");
    uProj_  = glGetUniformLocation(program_, "projection");
    uColor_ = glGetUniformLocation(program_, "objectColor");
}

void Renderer::initPrimitives() {
    float cubeVerts[] = {
        -1.f,-1.f,-1.f,  1.f, 1.f,-1.f,  1.f,-1.f,-1.f,
        -1.f,-1.f,-1.f, -1.f, 1.f,-1.f,  1.f, 1.f,-1.f,
        -1.f,-1.f, 1.f,  1.f,-1.f, 1.f,  1.f, 1.f, 1.f,
        -1.f,-1.f, 1.f,  1.f, 1.f, 1.f, -1.f, 1.f, 1.f,
        -1.f,-1.f,-1.f, -1.f,-1.f, 1.f, -1.f, 1.f, 1.f,
        -1.f,-1.f,-1.f, -1.f, 1.f, 1.f, -1.f, 1.f,-1.f,
         1.f,-1.f,-1.f,  1.f, 1.f, 1.f,  1.f,-1.f, 1.f,
         1.f,-1.f,-1.f,  1.f, 1.f,-1.f,  1.f, 1.f, 1.f,
        -1.f,-1.f,-1.f,  1.f,-1.f, 1.f,  1.f,-1.f,-1.f,
        -1.f,-1.f,-1.f, -1.f,-1.f, 1.f,  1.f,-1.f, 1.f,
        -1.f, 1.f,-1.f,  1.f, 1.f,-1.f,  1.f, 1.f, 1.f,
        -1.f, 1.f,-1.f,  1.f, 1.f, 1.f, -1.f, 1.f, 1.f
    };
    glGenVertexArrays(1, &cubeVAO_); 
    glGenBuffers(1, &cubeVBO_);
    glBindVertexArray(cubeVAO_);
    glBindBuffer(GL_ARRAY_BUFFER, cubeVBO_);
    glBufferData(GL_ARRAY_BUFFER, sizeof(cubeVerts), cubeVerts, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    cubeIndexCount_ = 36;
    
    glGenVertexArrays(1, &sphereVAO_); 
    glBindVertexArray(sphereVAO_);
    glGenBuffers(1, &sphereVBO_);
    glBindBuffer(GL_ARRAY_BUFFER, sphereVBO_);
    glBufferData(GL_ARRAY_BUFFER, 3 * sizeof(float), nullptr, GL_STATIC_DRAW); 
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    sphereIndexCount_ = 0; 
    glBindVertexArray(0); 
}

float Renderer::calculateRadius(double mass) const {
    double r_m = std::cbrt((3.0 * mass) / (4.0 * glm::pi<double>() * RHO_DEFAULT));
    return static_cast<float>(glm::max(r_m / SCALE_FACTOR, (double) MIN_RADIUS_SCALE));
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

    switch (mode_) {
    case Mode::Points: {
        glPointSize(8.0f);
        glBindVertexArray(sphereVAO_);
        for (const auto& obj : objs) {
            glm::vec3 pos = static_cast<glm::vec3>(obj.getPos()); 
            glm::vec4 col = obj.getColor(); 
            
            glm::mat4 M(1.0f);
            M = glm::translate(M, pos);
            glUniformMatrix4fv(uModel_, 1, GL_FALSE, glm::value_ptr(M));
            glUniform4f(uColor_, col.r, col.g, col.b, col.a);

            glDrawArrays(GL_POINTS, 0, 1);
        }
        glBindVertexArray(0);
        break;
    }

    case Mode::Cubes: {
        glBindVertexArray(cubeVAO_); 
        for (const auto& obj : objs) {
            float radius = calculateRadius(obj.getMass()); 
            glm::vec3 pos = static_cast<glm::vec3>(obj.getPos());
            glm::vec4 col = obj.getColor();
            glm::mat4 M(1.0f);
            M = glm::translate(M, pos);
            M = glm::scale(M, glm::vec3(radius)); 
            glUniformMatrix4fv(uModel_, 1, GL_FALSE, glm::value_ptr(M));
            glUniform4f(uColor_, col.r, col.g, col.b, col.a);
            glDrawArrays(GL_TRIANGLES, 0, cubeIndexCount_);
        }
        glBindVertexArray(0);
        break;
    }

    case Mode::Sphere:
    default: {
        glBindVertexArray(sphereVAO_); 
        for (const auto& obj : objs) {
            float radius = calculateRadius(obj.getMass()); 
            glm::vec3 pos = static_cast<glm::vec3>(obj.getPos());
            glm::vec4 col = obj.getColor();
            glm::mat4 M(1.0f);
            M = glm::translate(M, pos);
            M = glm::scale(M, glm::vec3(radius)); 
            glUniformMatrix4fv(uModel_, 1, GL_FALSE, glm::value_ptr(M));
            glUniform4f(uColor_, col.r, col.g, col.b, col.a);
            glDrawArrays(GL_TRIANGLES, 0, 36); 
        }
        glBindVertexArray(0);
        break;
    }
    }
}