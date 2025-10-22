#include "renderer.hpp"
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>
#include <cmath>

namespace {
    constexpr double kG = 6.6743e-11; // гравитац. постоянная для искажения сетки
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

void Renderer::setProjection(float fov_deg, float aspect, float znear, float zfar) {
    glUseProgram(program_);
    glm::mat4 P = glm::perspective(glm::radians(fov_deg), aspect, znear, zfar);
    glUniformMatrix4fv(uProj_, 1, GL_FALSE, glm::value_ptr(P));
}

void Renderer::updateView(const glm::vec3& pos,
                          const glm::vec3& front,
                          const glm::vec3& up) {
    glUseProgram(program_);
    glm::mat4 V = glm::lookAt(pos, pos + front, up);
    glUniformMatrix4fv(uView_, 1, GL_FALSE, glm::value_ptr(V));
}

void Renderer::updateGrid(float size, int divisions, const std::vector<Object>& objs) {
    auto verts = createGridVertices(size, divisions, objs);
    gridVertexCount_ = verts.size()/3;

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



void Renderer::drawObjects(const std::vector<Object>& objs) const {
    glUseProgram(program_);

    // ====== ОДНА ТОЧКА НА ОБЪЕКТ (АКТИВНО) ======
    // Рисуем каждый объект одной точкой. Никакой геометрии сферы/меша не используется.
    // Размер точки — настраиваемый:
    glPointSize(6.0f);

    for (const auto& obj : objs) {
        glm::mat4 M(1.0f);
        M = glm::translate(M, obj.position);
        glUniformMatrix4fv(uModel_, 1, GL_FALSE, glm::value_ptr(M));

        glUniform4f(uColor_, obj.color.r, obj.color.g, obj.color.b, obj.color.a);
        glBindVertexArray(obj.VAO);
        glDrawArrays(GL_POINTS, 0, 1);
    }
}




GLuint Renderer::compileProgram(const char* vs, const char* fs) {
    auto compile = [](GLenum type, const char* src){
        GLuint s = glCreateShader(type);
        glShaderSource(s, 1, &src, nullptr);
        glCompileShader(s);
        GLint ok=0; glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
        if(!ok){ char log[1024]; glGetShaderInfoLog(s,1024,nullptr,log);
            std::cerr << "Shader compile error: " << log << std::endl; }
        return s;
    };
    GLuint v = compile(GL_VERTEX_SHADER,   vs);
    GLuint f = compile(GL_FRAGMENT_SHADER, fs);

    GLuint p = glCreateProgram();
    glAttachShader(p, v); glAttachShader(p, f); glLinkProgram(p);
    GLint ok=0; glGetProgramiv(p, GL_LINK_STATUS, &ok);
    if(!ok){ char log[1024]; glGetProgramInfoLog(p,1024,nullptr,log);
        std::cerr << "Program link error: " << log << std::endl; }
    glDeleteShader(v); glDeleteShader(f);
    return p;
}

std::vector<float> Renderer::createGridVertices(float size, int divisions,
                                                const std::vector<Object>& objs) {
    std::vector<float> vertices;
    float step = size / divisions;
    float half = size * 0.5f;

    // X-параллельные
    for (int y=0; y<=divisions; ++y) {
        float yy = -half + y*step;
        for (int z=0; z<=divisions; ++z) {
            float zz = -half + z*step;
            for (int x=0; x<divisions; ++x) {
                float xs = -half + x*step;
                float xe = xs + step;
                vertices.insert(vertices.end(), {xs,yy,zz,  xe,yy,zz});
            }
        }
    }
    // Y-параллельные
    for (int x=0; x<=divisions; ++x) {
        float xx = -half + x*step;
        for (int z=0; z<=divisions; ++z) {
            float zz = -half + z*step;
            for (int y=0; y<divisions; ++y) {
                float ys = -half + y*step;
                float ye = ys + step;
                vertices.insert(vertices.end(), {xx,ys,zz,  xx,ye,zz});
            }
        }
    }
    // Z-параллельные
    for (int x=0; x<=divisions; ++x) {
        float xx = -half + x*step;
        for (int y=0; y<=divisions; ++y) {
            float yy = -half + y*step;
            for (int z=0; z<divisions; ++z) {
                float zs = -half + z*step;
                float ze = zs + step;
                vertices.insert(vertices.end(), {xx,yy,zs,  xx,yy,ze});
            }
        }
    }

    // Гравитационное смещение
    for (size_t i=0; i<vertices.size(); i+=3) {
        glm::vec3 p(vertices[i], vertices[i+1], vertices[i+2]);
        glm::vec3 disp(0.0f);
        for (const auto& o : objs) {
            glm::vec3 d = o.GetPos() - p;
            float r = glm::length(d);
            float r_m = r * 1000.0f;
            if (r_m < 1e-5f) continue; // защита от деления на ноль
            float strength = float((kG * o.mass) / double(r_m*r_m));
            glm::vec3 dir = (r > 0.0f) ? d / r : glm::vec3(0);
            glm::vec3 one = dir * strength;
            if (! (r + glm::length(one) < o.radius * 1000.0f)) {
                disp += one * (2.0f / glm::max(r, 1e-4f));
            }
        }
        p += disp;
        vertices[i] = p.x; vertices[i+1] = p.y; vertices[i+2] = p.z;
    }
    return vertices;
}
