#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <iostream>
#include <optional>
#include <cmath>
#include <algorithm>

#include "object.hpp"
#include "renderer.hpp"
#include "control.hpp"

bool running = true;
bool pause   = false;

glm::vec3 cameraPos   = glm::vec3(0.0f, 0.0f, 1.0f);
glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
glm::vec3 cameraUp    = glm::vec3(0.0f, 1.0f, 0.0f);

float lastX = 400.0f, lastY = 300.0f;
float yaw   = -90.0f;
float pitch = 0.0f;

float deltaTime = 0.0f;
float lastFrame = 0.0f;

float timeScale = 1.0f; // переменная для ускорения/замедления времени

float grid_size2  = 400.0f;
int   vert_count2 = 10;

const double G = 6.6743e-11;
float initMass = 5.0f * std::pow(10.0f, 20.0f) / 5.0f;
char title[128];

std::vector<Object> objs = {};

// Шейдеры
const char* vertexShaderSource = R"glsl(
#version 330 core
layout(location=0) in vec3 aPos;
uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
void main() {
    gl_Position = projection * view * model * vec4(aPos, 1.0);
}
)glsl";

const char* fragmentShaderSource = R"glsl(
#version 330 core
out vec4 FragColor;
uniform vec4 objectColor;
void main() {
    FragColor = objectColor;
}
)glsl";

GLFWwindow* StartGLU();

int main() {

    GLFWwindow* window = StartGLU();
    if (!window) {
        std::cerr << "Window or OpenGL context creation failed.\n";
        return -1;
    }

    Renderer renderer(800, 600, vertexShaderSource, fragmentShaderSource);
    renderer.setProjection(65.0f, 800.0f/600.0f, 0.1f, 1200.0f);

    cameraPos = glm::vec3(0.0f, 50.0f, 250.0f);


objs = {
    Object(glm::vec3(0.0f, 0.0f, 0.0f),
           glm::vec3(0.0f, 0.0f, 0.0f),
           initMass * 1200.0, 141000.0f, std::nullopt),

    Object(glm::vec3( 60.0f, 0.0f,   0.0f),
           glm::vec3(  100.0f, 0.0f, 0.0f),
           initMass * 1.0, 1410.0f, std::nullopt),
    Object(glm::vec3( -120.0f, 0.0f,   0.0f),
           glm::vec3(  0.0f, 100.0f, 0.0f),
           initMass * 1.0, 1410.0f, std::nullopt),
    Object(glm::vec3( -180.0f, 0.0f,   0.0f),
           glm::vec3(  0.0f, 100.0f, 0.0f),
           initMass * 1.0, 1410.0f, std::nullopt)
};


    // Управление 
    Control control(window, objs,
                    cameraPos, cameraFront, cameraUp,
                    deltaTime, timeScale, pause, running,
                    yaw, pitch, lastX, lastY,
                    initMass);
    control.attach();

    while (!glfwWindowShouldClose(window) && running) {
        float currentFrame = glfwGetTime();
        float dtReal = currentFrame - lastFrame;
        deltaTime = dtReal * timeScale;
        lastFrame = currentFrame;

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        renderer.updateView(cameraPos, cameraFront, cameraUp);

        // Увел. массы объекта при зажатой клавиши
        if (!objs.empty() && objs.back().Initalizing) {
            if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS) {
                objs.back().mass *= 1.0f + 1.0f * deltaTime;
                objs.back().radius = std::pow(
                    (3 * objs.back().mass / objs.back().density) / (4 * 3.14159265359f),
                    1.0f/3.0f
                ) / 100000.0f;
                objs.back().UpdateVertices();
            }
        }

        for (size_t i = 0; i < objs.size(); ++i) {
            Object& obj = objs[i];
            if (obj.Initalizing) {
                obj.radius = std::pow((3 * obj.mass / obj.density) / (4 * 3.14159265359f), 1.0f/3.0f) / 100000.0f;
                obj.UpdateVertices();
            }

            for (size_t j = i + 1; j < objs.size(); ++j) {
                Object& obj2 = objs[j];
                if (obj.Initalizing || obj2.Initalizing) {
                    continue;
                }

                glm::vec3 delta = obj2.GetPos() - obj.GetPos();
                float distance = std::sqrt(delta.x * delta.x + delta.y * delta.y + delta.z * delta.z);
                if (distance <= 0.0f) {
                    continue;
                }

                glm::vec3 dir = delta / distance;
                float combinedRadius = obj.radius + obj2.radius;

                float effectiveDistance = std::max(distance, combinedRadius);
                float dist_m = effectiveDistance * 1000.0f;
                double F = (G * obj.mass * obj2.mass) / (dist_m * dist_m);
                float acc1 = static_cast<float>(F / obj.mass);
                float acc2 = static_cast<float>(F / obj2.mass);
                glm::vec3 accObj = dir * acc1;
                glm::vec3 accObj2 = -dir * acc2;

                if (!pause) {
                    obj.accelerate(accObj.x, accObj.y, accObj.z, deltaTime);
                    obj2.accelerate(accObj2.x, accObj2.y, accObj2.z, deltaTime);
                }

                if (distance < combinedRadius) {
                    glm::vec3 normal = dir;
                    glm::vec3 relativeVelocity = obj.velocity - obj2.velocity;
                    float relVelAlongNormal = glm::dot(relativeVelocity, normal);

                    if (relVelAlongNormal < 0.0f) {
                        float restitution = 0.8f;
                        float invMass1 = 1.0f / obj.mass;
                        float invMass2 = 1.0f / obj2.mass;
                        float impulseScalar = -(1.0f + restitution) * relVelAlongNormal / (invMass1 + invMass2);
                        glm::vec3 impulse = impulseScalar * normal;
                        obj.velocity += impulse * invMass1;
                        obj2.velocity -= impulse * invMass2;
                    }

                    float penetration = combinedRadius - distance;
                    if (penetration > 0.0f) {
                        float invMass1 = 1.0f / obj.mass;
                        float invMass2 = 1.0f / obj2.mass;
                        float invMassSum = invMass1 + invMass2;
                        if (invMassSum > 0.0f) {
                            glm::vec3 correction = normal * (penetration / invMassSum);
                            obj.position -= correction * invMass1;
                            obj2.position += correction * invMass2;
                        }
                    }
                }
            }

            if (!pause) {
                obj.UpdatePos(deltaTime);
            }
        }

        // Сетка
        renderer.updateGrid(grid_size2, vert_count2, objs);
        renderer.drawGrid();

        // Отрисовка 
        renderer.drawObjects(objs);

        float fps = (dtReal > 0.0f) ? 1.0f / dtReal : 0.0f;
        std::snprintf(title, sizeof(title),
              "3D_TEST | timeScale: %.2fx | FPS: %.0f | Objects =: %zu",
              timeScale, fps, objs.size());
        glfwSetWindowTitle(window, title);
        glfwSwapBuffers(window);

        glfwPollEvents();
    }

    for (auto& obj : objs) {
        glDeleteVertexArrays(1, &obj.VAO);
        glDeleteBuffers(1, &obj.VBO);
    }
    glfwTerminate();
    return 0;
}

// Инициализируем контекст openGL
GLFWwindow* StartGLU() {

    glfwInitHint(GLFW_PLATFORM, GLFW_ANY_PLATFORM);

    if (!glfwInit()) {
        std::cout << "Failed to initialize GLFW\n";
        return nullptr;
    }

    GLFWwindow* window = glfwCreateWindow(800, 600, "3D_TEST", nullptr, nullptr);
    if (!window) {
        std::cerr << "Failed to create GLFW window.\n";
        glfwTerminate();
        return nullptr;
    }
    glfwMakeContextCurrent(window);

    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW.\n";
        glfwTerminate();
        return nullptr;
    }

    glEnable(GL_DEPTH_TEST);
    glViewport(0, 0, 800, 600);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    return window;
}