#pragma once
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <vector>
#include "object.hpp"

enum class RenderMode { Sphere, Points, Cubes };

class Renderer {
public:
    Renderer(int w, int h, const char* vs, const char* fs);
    ~Renderer();

    void setProjection(float fov_deg, float aspect, float znear, float zfar);
    void updateView(const glm::vec3& pos, const glm::vec3& front, const glm::vec3& up);
    void drawObjects(const std::vector<Object>& objs) const;
    
    void setRenderMode(RenderMode mode) { mode_ = mode; }

private:
    GLuint compileProgram(const char* vs, const char* fs);
    
    // Инициализаторы геометрии
    void initSphereGeometry();
    void initCubeGeometry();  
    void initPointGeometry(); 

    GLuint program_;
    GLint uModel_, uView_, uProj_, uColor_;
    
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
};