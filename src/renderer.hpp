#pragma once
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <vector>
#include "object.hpp"

class Renderer {
public:
    Renderer(int width, int height,
             const char* vertexSrc,
             const char* fragmentSrc);
    ~Renderer();

    void setProjection(float fov_deg, float aspect, float znear, float zfar);
    void updateView(const glm::vec3& pos,
                    const glm::vec3& front,
                    const glm::vec3& up);

    void updateGrid(float size, int divisions, const std::vector<Object>& objs);

    void drawGrid() const;
    void drawObjects(const std::vector<Object>& objs) const;

    static float ComputeFar(glm::vec3& camPos, std::vector<Object>& objs);

private:
    GLuint program_ = 0;
    GLint  uModel_ = -1, uView_ = -1, uProj_ = -1, uColor_ = -1;

    GLuint gridVAO_ = 0, gridVBO_ = 0;
    size_t gridVertexCount_ = 0;

    static GLuint compileProgram(const char* vs, const char* fs);
    static std::vector<float> createGridVertices(float size, int divisions,
                                                 const std::vector<Object>& objs);
};
