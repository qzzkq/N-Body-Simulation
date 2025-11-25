#include "control.hpp"
#include "constants.hpp"
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
}

// Статические обертки (без изменений)
void Control::KeyCB(GLFWwindow* w, int key, int scancode, int action, int mods) {
    static_cast<Control*>(glfwGetWindowUserPointer(w))->onKey(key, scancode, action, mods);
}
void Control::MouseButtonCB(GLFWwindow* w, int button, int action, int mods) {
    static_cast<Control*>(glfwGetWindowUserPointer(w))->onMouseButton(button, action, mods);
}
void Control::ScrollCB(GLFWwindow* w, double xoffset, double yoffset) {
    static_cast<Control*>(glfwGetWindowUserPointer(w))->onScroll(xoffset, yoffset);
}
void Control::CursorPosCB(GLFWwindow* w, double xpos, double ypos) {
    static_cast<Control*>(glfwGetWindowUserPointer(w))->onCursorPos(xpos, ypos);
}

// Реальная логика клавиатуры
void Control::onKey(int key, int /*scancode*/, int action, int /*mods*/) {
    if (action != GLFW_PRESS && action != GLFW_REPEAT) return;

    float cameraSpeed = 1000.0f * deltaTime_;
    if (glfwGetKey(window_, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) cameraSpeed *= 10.0f;
    
    // Камера (используются ссылки на AppState)
    if (key == GLFW_KEY_W) cameraPos_ += cameraSpeed * cameraFront_;
    if (key == GLFW_KEY_S) cameraPos_ -= cameraSpeed * cameraFront_;
    if (key == GLFW_KEY_A) cameraPos_ -= glm::normalize(glm::cross(cameraFront_, cameraUp_)) * cameraSpeed;
    if (key == GLFW_KEY_D) cameraPos_ += glm::normalize(glm::cross(cameraFront_, cameraUp_)) * cameraSpeed;
    if (key == GLFW_KEY_SPACE) cameraPos_ += cameraUp_ * cameraSpeed;
    if (key == GLFW_KEY_LEFT_CONTROL) cameraPos_ -= cameraUp_ * cameraSpeed;

    // Управление приложением (используются ссылки на AppState)
    if (key == GLFW_KEY_ESCAPE) running_ = false;
    if (key == GLFW_KEY_P) pause_ = !pause_;
    if (key == GLFW_KEY_UP) timeScale_ = std::min(timeScale_ + 0.1f, 5.0f);
    if (key == GLFW_KEY_DOWN) timeScale_ = std::max(timeScale_ - 0.1f, 0.1f);
}

// Реальная логика мыши (создание объекта)
void Control::onMouseButton(int button, int action, int /*mods*/) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) {
            // Создание нового объекта
            objs_.emplace_back(
                glm::dvec3(cameraPos_), // Pos
                glm::dvec3(0.0),        // Vel
                initMass_,              // Mass
                glm::dvec3(0.0)         // Acc
            );
            objs_.back().Initalizing = true;
        }
        if (action == GLFW_RELEASE) {
            objs_.back().Initalizing = false;
            objs_.back().Launched = true;
            
            // Расчет начальной скорости на основе смещения
            // (Используем геттеры/сеттеры!)
            double x_cursor = lastX_;
            double y_cursor = lastY_;
            // Нужно получить координаты для расчета вектора скорости
            // В твоем старом коде здесь была логика, которую сложно восстановить без полных данных.
            // Предположим, что скорость задается на основе позиции камеры
            // glm::dvec3 newVel = (desired_pos - objs_.back().getPos()) * SCALE;
            // objs_.back().setVel(newVel);
        }
    }
}

void Control::onScroll(double /*xoffset*/, double yoffset) {
    float cameraSpeed = 5000.0f * deltaTime_;
    if (yoffset > 0) cameraPos_ += cameraSpeed * cameraFront_;
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

    if (pitch_ > 89.0f) pitch_ = 89.0f;
    if (pitch_ < -89.0f) pitch_ = -89.0f;

    // Обновление cameraFront (используется ссылка на AppState)
    glm::vec3 front;
    front.x = cos(glm::radians(yaw_)) * cos(glm::radians(pitch_));
    front.y = sin(glm::radians(pitch_));
    front.z = sin(glm::radians(yaw_)) * cos(glm::radians(pitch_));
    cameraFront_ = glm::normalize(front);
}
