#include "windowHandler.hpp"

WindowHandler::WindowHandler() {
    this->window = nullptr; 
}

double WindowHandler::getTime() {
    return glfwGetTime(); 
}

int WindowHandler::createWindow() {
    if (!glfwInit()) {
        return 0;
    }
    glfwInitHint(GLFW_PLATFORM, GLFW_PLATFORM_X11);
    GLFWwindow* window = glfwCreateWindow(800, 600, "3D_TEST", nullptr, nullptr);
    if (!window) {
        glfwTerminate();
        return 0;
    }

    glfwMakeContextCurrent(window);

    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK) {
        glfwTerminate();
        return 0;
    }

    glEnable(GL_DEPTH_TEST);
    glViewport(0, 0, 800, 600);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    this->window = window;
    return 1; 
}

void WindowHandler::setTitle(char* str) {
    glfwSetWindowTitle(this->window, str);
}

void WindowHandler::cleanResources() {
    glfwTerminate();
}

void WindowHandler::show() {
    glfwSwapBuffers(this->window);
    glfwPollEvents();
}

GLFWwindow* WindowHandler::getWindowPointer() {
    return this->window; 
}

int WindowHandler::closeFlag() {
    return glfwWindowShouldClose(this->window); 
}

void WindowHandler::cleanWindow() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}