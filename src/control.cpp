#include "control.hpp"
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <cmath>

Control::Control(GLFWwindow* window,
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
                 float& initMass)
    : window_(window),
      objs_(objs),
      cameraPos_(cameraPos),
      cameraFront_(cameraFront),
      cameraUp_(cameraUp),
      deltaTime_(deltaTime),
      timeScale_(timeScale),
      pause_(pause),
      running_(running),
      yaw_(yaw),
      pitch_(pitch),
      lastX_(lastX),
      lastY_(lastY),
      initMass_(initMass)
{}

void Control::attach() {
    glfwSetWindowUserPointer(window_, this);
    glfwSetKeyCallback(window_, &Control::KeyCB);
    glfwSetMouseButtonCallback(window_, &Control::MouseButtonCB);
    glfwSetScrollCallback(window_, &Control::ScrollCB);
    glfwSetCursorPosCallback(window_, &Control::CursorPosCB);
    glfwSetInputMode(window_, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
}

/* static */ void Control::KeyCB(GLFWwindow* w, int key, int sc, int action, int mods) {
    auto* self = static_cast<Control*>(glfwGetWindowUserPointer(w));
    if (self) self->onKey(key, sc, action, mods);
}
/* static */ void Control::MouseButtonCB(GLFWwindow* w, int button, int action, int mods) {
    auto* self = static_cast<Control*>(glfwGetWindowUserPointer(w));
    if (self) self->onMouseButton(button, action, mods);
}
/* static */ void Control::ScrollCB(GLFWwindow* w, double xo, double yo) {
    auto* self = static_cast<Control*>(glfwGetWindowUserPointer(w));
    if (self) self->onScroll(xo, yo);
}
/* static */ void Control::CursorPosCB(GLFWwindow* w, double x, double y) {
    auto* self = static_cast<Control*>(glfwGetWindowUserPointer(w));
    if (self) self->onCursorPos(x, y);
}

void Control::onKey(int key, int /*scancode*/, int action, int mods) {
    float cameraSpeed = 1000.0f * deltaTime_;
    bool shiftPressed = (mods & GLFW_MOD_SHIFT) != 0;

    // Камера WASD + Space/Shift
    if (glfwGetKey(window_, GLFW_KEY_W) == GLFW_PRESS) cameraPos_ += cameraSpeed * cameraFront_;
    if (glfwGetKey(window_, GLFW_KEY_S) == GLFW_PRESS) cameraPos_ -= cameraSpeed * cameraFront_;
    if (glfwGetKey(window_, GLFW_KEY_A) == GLFW_PRESS) cameraPos_ -= cameraSpeed * glm::normalize(glm::cross(cameraFront_, cameraUp_));
    if (glfwGetKey(window_, GLFW_KEY_D) == GLFW_PRESS) cameraPos_ += cameraSpeed * glm::normalize(glm::cross(cameraFront_, cameraUp_));
    if (glfwGetKey(window_, GLFW_KEY_SPACE) == GLFW_PRESS)      cameraPos_ += cameraSpeed * cameraUp_;
    if (glfwGetKey(window_, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) cameraPos_ -= cameraSpeed * cameraUp_;

    // Ускорение/замедление времени
    if (glfwGetKey(window_, GLFW_KEY_EQUAL)   == GLFW_PRESS) timeScale_ = std::min(timeScale_ * 1.1f, 100.0f);
    if (glfwGetKey(window_, GLFW_KEY_MINUS)   == GLFW_PRESS) timeScale_ = std::max(timeScale_ / 1.1f, 0.05f);
    
    // Пауза на K
    if (glfwGetKey(window_, GLFW_KEY_K) == GLFW_PRESS)   pause_ = true;
    if (glfwGetKey(window_, GLFW_KEY_K) == GLFW_RELEASE) pause_ = false;

    // Выход на Q
    if (glfwGetKey(window_, GLFW_KEY_Q) == GLFW_PRESS) {
        glfwTerminate();
        glfwSetWindowShouldClose(window_, GLFW_TRUE);
        running_ = false;
    }

    // Движение создаваемого объекта стрелками
    if (!objs_.empty() && objs_.back().Initalizing) {
        if (key == GLFW_KEY_UP && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
            if (!shiftPressed) objs_.back().position.y += 0.5f;
            objs_.back().position.z += 0.5f;
        }
        if (key == GLFW_KEY_DOWN && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
            if (!shiftPressed) objs_.back().position.y -= 0.5f;
            objs_.back().position.z -= 0.5f;
        }
        if (key == GLFW_KEY_RIGHT && (action == GLFW_PRESS || action == GLFW_REPEAT))
            objs_.back().position.x += 0.5f;
        if (key == GLFW_KEY_LEFT  && (action == GLFW_PRESS || action == GLFW_REPEAT))
            objs_.back().position.x -= 0.5f;
    }
}

void Control::onMouseButton(int button, int action, int /*mods*/) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) {
            objs_.emplace_back(glm::vec3(0.0f),
                               glm::vec3(0.0f),
                               initMass_, 1000.0f, std::nullopt);
            objs_.back().Initalizing = true;
        }
        if (action == GLFW_RELEASE) {
            objs_.back().Initalizing = false;
            objs_.back().Launched = true;
        }
    }
}

void Control::onScroll(double /*xoffset*/, double yoffset) {
    float cameraSpeed = 5000.0f * deltaTime_;
    if (yoffset > 0)      cameraPos_ += cameraSpeed * cameraFront_;
    else if (yoffset < 0) cameraPos_ -= cameraSpeed * cameraFront_;
}

void Control::onCursorPos(double xpos, double ypos) {
    float xoffset = static_cast<float>(xpos - lastX_);
    float yoffset = static_cast<float>(lastY_ - ypos);
    lastX_ = static_cast<float>(xpos);
    lastY_ = static_cast<float>(ypos);

    float sensitivity = 0.1f;
    xoffset *= sensitivity;
    yoffset *= sensitivity;

    yaw_   += xoffset;
    pitch_ += yoffset;
    if (pitch_ > 89.0f)  pitch_ = 89.0f;
    if (pitch_ < -89.0f) pitch_ = -89.0f;

    glm::vec3 front;
    front.x = cos(glm::radians(yaw_)) * cos(glm::radians(pitch_));
    front.y = sin(glm::radians(pitch_));
    front.z = sin(glm::radians(yaw_)) * cos(glm::radians(pitch_));
    cameraFront_ = glm::normalize(front);
}
