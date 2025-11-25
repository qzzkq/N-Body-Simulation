#pragma once
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <vector>
#include "object.hpp"

class Renderer {
public:
    // режим визуализации, выбираем его из main
    enum class Mode {Sphere, Points, Cubes};

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

    GLuint sphereVAO_ = 0, sphereVBO_ = 0; 
    GLuint cubeVAO_ = 0, cubeVBO_ = 0;
    size_t cubeIndexCount_ = 0;
    size_t sphereIndexCount_ = 0;

    // выбранный режим отрисовки
    Mode mode_ = Mode::Sphere;

    float calculateRadius(double mass) const; // подсчёт размера из массы 
    compileProgram(const char*, const char*); 
    void initShader();
    void initPrimitives();  

};
