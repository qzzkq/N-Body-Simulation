#pragma once
struct GLFWwindow;
#include <glm/glm.hpp>
#include <vector>
#include "object.hpp"

// Управляет вводом, камерой и созданием объектов.
class Control {
public:
    Control(GLFWwindow* window,
            std::vector<Object>& objs,
            glm::vec3& cameraPos,
            glm::vec3& cameraFront,
            glm::vec3& cameraUp,
            float& deltaTime,
            float& timeScale,
            bool& pause,
            bool& running,
            float& yaw,
            float& pitch,
            float& lastX,
            float& lastY,
            float& initMass);

    void attach();

private:
    // Статические обёртки для GLFW
    static void KeyCB(GLFWwindow* w, int key, int scancode, int action, int mods);
    static void MouseButtonCB(GLFWwindow* w, int button, int action, int mods);
    static void ScrollCB(GLFWwindow* w, double xoffset, double yoffset);
    static void CursorPosCB(GLFWwindow* w, double xpos, double ypos);

    // Реальная логика обработчиков
    void onKey(int key, int scancode, int action, int mods);
    void onMouseButton(int button, int action, int mods);
    void onScroll(double xoffset, double yoffset);
    void onCursorPos(double xpos, double ypos);

    // Ссылки на состояние
    GLFWwindow* window_;
    std::vector<Object>& objs_;
    glm::vec3& cameraPos_;
    glm::vec3& cameraFront_;
    glm::vec3& cameraUp_;
    float& deltaTime_;
    float& timeScale_;
    bool& pause_;
    bool& running_;
    float& yaw_;
    float& pitch_;
    float& lastX_;
    float& lastY_;
    float& initMass_;
};
