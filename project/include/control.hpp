#pragma once
struct GLFWwindow;
#include <glm/glm.hpp>
#include <vector>
#include "object.hpp"
#include "camera.hpp"
#include "state.hpp"

// Управляет вводом, камерой и созданием объектов.
class Control {
public:
    Control(GLFWwindow* window, 
            std::vector<Object>& objs, 
            Camera& cam, 
            SimState& state);

    void attach();

private:
    // Статические обёртки для GLFW
    static void KeyCB(GLFWwindow* w, int key, int scancode, int action, int mods);
    static void ScrollCB(GLFWwindow* w, double xoffset, double yoffset);
    static void CursorPosCB(GLFWwindow* w, double xpos, double ypos);

    // Реальная логика обработчиков
    void onKey(int key, int scancode, int action, int mods);
    void onScroll(double xoffset, double yoffset);
    void onCursorPos(double xpos, double ypos);

    // Ссылки на состояние
    GLFWwindow* window_;
    std::vector<Object>& objs_;
    Camera& camera_;
    SimState& state_;
};
