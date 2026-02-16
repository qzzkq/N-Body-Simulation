#ifndef CAMERAINC
#define CAMERAINC

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

typedef struct Camera {
    glm::vec3 pos   = glm::vec3(0.0f, 0.0f, 1.0f);
    glm::vec3 front = glm::vec3(0.0f, 0.0f, -1.0f);
    glm::vec3 up    = glm::vec3(0.0f, 1.0f, 0.0f);
    float yaw   = -90.0f;
    float pitch = 0.0f;
    float lastX = 400.0f;
    float lastY = 300.0f;
    float fov   = 45.0f;

    glm::mat4 getViewMatrix() const {
        return glm::lookAt(pos, pos + front, up);
    }
} Camera;

#endif 