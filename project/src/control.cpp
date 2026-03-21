#include "control.hpp"
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <algorithm>
#include <cmath>

Control::Control(GLFWwindow* window, 
                std::vector<Object>& objs, 
                Camera& cam, 
                SimState& state)
    : window_(window),
      objs_(objs),
      camera_(cam), 
      state_(state) 
{}

void Control::attach() {
    glfwSetWindowUserPointer(window_, this);
    glfwSetKeyCallback(window_, &Control::KeyCB);
    glfwSetScrollCallback(window_, &Control::ScrollCB);
    glfwSetCursorPosCallback(window_, &Control::CursorPosCB);
    glfwSetInputMode(window_, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
}

/* static */ void Control::KeyCB(GLFWwindow* w, int key, int sc, int action, int mods) {
    auto* self = static_cast<Control*>(glfwGetWindowUserPointer(w));
    if (self) self->onKey(key, sc, action, mods);
}

/* static */ void Control::ScrollCB(GLFWwindow* w, double xo, double yo) {
    auto* self = static_cast<Control*>(glfwGetWindowUserPointer(w));
    if (self) self->onScroll(xo, yo);
}
/* static */ void Control::CursorPosCB(GLFWwindow* w, double x, double y) {
    auto* self = static_cast<Control*>(glfwGetWindowUserPointer(w));
    if (self) self->onCursorPos(x, y);
}

void Control::updateCameraFromKeys() {
    float cameraSpeed = 25.0f * cameraMoveScale_ * state_.deltaTime;
    if (cameraSpeed == 0.0f) return;

    // Камера WASD + Space / Ctrl (вниз). Shift зарезервирован под «заморозку времени» в main/replay.
    if (glfwGetKey(window_, GLFW_KEY_W) == GLFW_PRESS) camera_.pos += cameraSpeed * camera_.front;
    if (glfwGetKey(window_, GLFW_KEY_S) == GLFW_PRESS) camera_.pos -= cameraSpeed * camera_.front;
    if (glfwGetKey(window_, GLFW_KEY_A) == GLFW_PRESS) camera_.pos -= cameraSpeed * glm::normalize(glm::cross(camera_.front, camera_.up));
    if (glfwGetKey(window_, GLFW_KEY_D) == GLFW_PRESS) camera_.pos += cameraSpeed * glm::normalize(glm::cross(camera_.front, camera_.up));
    if (glfwGetKey(window_, GLFW_KEY_SPACE) == GLFW_PRESS) camera_.pos += cameraSpeed * camera_.up;
    if (glfwGetKey(window_, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS) camera_.pos -= cameraSpeed * camera_.up;
}

void Control::onKey(int key, int /*scancode*/, int action, int mods) {
    const bool shift = (mods & GLFW_MOD_SHIFT) != 0;
    const bool repeat = (action == GLFW_PRESS || action == GLFW_REPEAT);

    // Shift + − / Shift + =: скорость перемещения камеры (WASD, колёсико). Без Shift — время симуляции.
    if (repeat && shift) {
        if (key == GLFW_KEY_MINUS || key == GLFW_KEY_KP_SUBTRACT) {
            cameraMoveScale_ = std::max(cameraMoveScale_ / 1.1f, 0.01f);
            return;
        }
        if (key == GLFW_KEY_EQUAL || key == GLFW_KEY_KP_ADD) {
            cameraMoveScale_ *= 1.1f;
            return;
        }
    }

    if (action == GLFW_PRESS && !shift) {
        if (key == GLFW_KEY_EQUAL || key == GLFW_KEY_KP_ADD) {
            state_.timeScale *= 1.1f;
        } else if (key == GLFW_KEY_MINUS || key == GLFW_KEY_KP_SUBTRACT) {
            state_.timeScale = std::max(state_.timeScale / 1.1f, 0.05f);
        }
    }

    // Пауза на K
    if (key == GLFW_KEY_K) {
        if (action == GLFW_PRESS) state_.pause = true;
        if (action == GLFW_RELEASE) state_.pause = false;
    }

    // Выход на Q
    if (key == GLFW_KEY_Q && action == GLFW_PRESS) {
        glfwTerminate();
        glfwSetWindowShouldClose(window_, GLFW_TRUE);
        state_.running = false;
    }

    //регулировка размера окна 
    // уменьшение
    if (key == GLFW_KEY_LEFT_BRACKET && action == GLFW_PRESS) {
        int w = 0, h = 0;  
        glfwGetWindowSize(window_, &w, &h);
        w = (int) (w * 0.9);
        h = (int) (h * 0.9); 
        glfwSetWindowSize(window_, w, h); // установление логического размера 
        int frW = 0, frH = 0; 
        glfwGetFramebufferSize(window_, &frW, &frH); 
        glViewport(0, 0, frW, frH); 
    }

    //увеличение
    if (key == GLFW_KEY_RIGHT_BRACKET && action == GLFW_PRESS) {
        int w = 0; 
        int h = 0; 
        glfwGetWindowSize(window_, &w, &h);
        w = (int) (w * 1.1);
        h = (int) (h * 1.1); 
        glfwSetWindowSize(window_, w, h); // установление логического размера
        int frW = 0, frH = 0; 
        glfwGetFramebufferSize(window_, &frW, &frH); 
        glViewport(0, 0, frW, frH); 
    }

}

void Control::onScroll(double /*xoffset*/, double yoffset) {
    float cameraSpeed = 100.0f * cameraMoveScale_ * state_.deltaTime;
    if (yoffset > 0)      camera_.pos += cameraSpeed *  camera_.front;
    else if (yoffset < 0) camera_.pos -= cameraSpeed *  camera_.front;
}

void Control::onCursorPos(double xpos, double ypos) {
    float xoffset = static_cast<float>(xpos - camera_.lastX);
    float yoffset = static_cast<float>(camera_.lastY - ypos);
    camera_.lastX = static_cast<float>(xpos);
    camera_.lastY = static_cast<float>(ypos);

    float sensitivity = 0.1f;
    xoffset *= sensitivity;
    yoffset *= sensitivity;

    camera_.yaw   += xoffset;
    camera_.pitch += yoffset;
    if (camera_.pitch > 89.0f)  camera_.pitch = 89.0f;
    if (camera_.pitch < -89.0f) camera_.pitch = -89.0f;

    glm::vec3 front;
    front.x = cos(glm::radians(camera_.yaw)) * cos(glm::radians(camera_.pitch));
    front.y = sin(glm::radians(camera_.pitch));
    front.z = sin(glm::radians(camera_.yaw)) * cos(glm::radians(camera_.pitch));
    camera_.front = glm::normalize(front);
}
