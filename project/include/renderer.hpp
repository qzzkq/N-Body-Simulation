#pragma once
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <GLFW/glfw3.h>
#include <vector>
#include "object.hpp"
#include "camera.hpp" 

enum class RenderMode { Sphere, Points, Cubes };

class Renderer {
public:
    Renderer();
    ~Renderer();

    bool init(int width,  int height, const char* title, bool fullscreen, bool maximized); 

    void setProjection(float fov_deg, float znear, float zfar);
    void updateView(const Camera& camera);
    void drawObjects(const std::vector<Object>& objs) const;
    
    void setRenderMode(RenderMode mode) { mode_ = mode; }

    GLFWwindow* getWindow() const { return window_; }
    
    void renderFrame(const std::vector<Object>& objs, const Camera& cam);
    
    void resizeWindow(int w, int h);

    void drawTrails(const std::vector<Object>& objs) const; 

private:
    GLuint compileProgram(const char* vs, const char* fs);
    
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
    GLint uModel_, uView_, uProj_, uColor_;
    
    GLuint trailProgram_;

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

    //VAO и VBO для трэйлов
    GLuint trailVAO_ = 0;
    GLuint trailVBO_ = 0; 

    glm::mat4 viewMatrix_; 
    glm::mat4 projectionMatrix_; 
    
    bool successInit_ = false; 
};