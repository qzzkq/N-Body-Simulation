#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <iostream>
#include <vector>
#include <string>
#include <cstdio>
#include <cmath>

#include "object.hpp"
#include "renderer.hpp"
#include "data.hpp"
#include "H5Cpp.h"

// ==== те же шейдеры, что и в main.cpp ====

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

// ==== Копия StartGLU из main.cpp ====

GLFWwindow* StartGLU() {
    if (!glfwInit()) {
        std::cout << "Failed to initialize GLFW\n";
        return nullptr;
    }
    glfwInitHint(GLFW_PLATFORM, GLFW_PLATFORM_X11);
    GLFWwindow* window = glfwCreateWindow(800, 600, "Replay", nullptr, nullptr);
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

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Использование: " << argv[0] << " output.h5\n";
        return 1;
    }

    std::string fileName = argv[1];

    // --- DETECT: считаем количество фреймов по наличию датасетов frame_XXXXXX ---

    std::size_t numFrames = 0;
    std::size_t numBodies = 0;
    std::vector<Particle> firstFrameParts;

    while (true) {
        char dsetName[64];
        std::snprintf(dsetName, sizeof(dsetName), "frame_%06zu", numFrames);

        try {
            auto parts = Reader(fileName, dsetName);
            if (parts.empty()) {
                break;
            }

            if (numFrames == 0) {
                numBodies = parts.size();
                firstFrameParts = parts; // копируем начальный кадр
            } else {
                if (parts.size() != numBodies) {
                    std::cerr << "Кадр " << dsetName
                              << " имеет " << parts.size()
                              << " тел, а первый кадр имел " << numBodies << "\n";
                    return 1;
                }
            }

            ++numFrames;
        }
        catch (const H5::Exception&) {
            // датасета с таким именем нет → все фреймы закончились
            break;
        }
        catch (const std::exception& e) {
            std::cerr << "Ошибка при чтении " << dsetName << ": " << e.what() << "\n";
            return 1;
        }
    }

    if (numFrames == 0 || numBodies == 0) {
        std::cerr << "В файле " << fileName
                  << " не найдено ни одного кадра frame_XXXXXX.\n";
        return 1;
    }

    std::cout << "Файл: " << fileName
              << "\n  num_bodies (детект) = " << numBodies
              << "\n  num_frames (детект) = " << numFrames << "\n";

    // --- OpenGL / окно ---
    GLFWwindow* window = StartGLU();
    if (!window) {
        return 1;
    }

    Renderer renderer(800, 600, vertexShaderSource, fragmentShaderSource);
    renderer.setProjection(65.0f, 800.0f / 600.0f, 8.3f, 100000.0f);

    glm::vec3 cameraPos   = glm::vec3(0.0f, 50.0f, 250.0f);
    glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
    glm::vec3 cameraUp    = glm::vec3(0.0f, 1.0f,  0.0f);

    // --- Создаём объекты по первому кадру ---

    std::vector<Object> objs;
    objs.reserve(numBodies);
    for (const auto& p : firstFrameParts) {
        Object o(p.position, p.velocity,
                 p.mass,            // mass
                 1410.0f,           // плотность как в spawnSystem
                 std::nullopt);     // radius потом зададим

        o.Initalizing = false;
        if (std::isfinite(p.radius) && p.radius > 0.0) {
            o.radius = static_cast<float>(p.radius);
        } else {
            o.radius = 1.0f;
        }
        o.UpdateVertices();
        objs.push_back(std::move(o));
    }

    std::size_t currentFrame = 0;
    double lastSwitch = glfwGetTime();
    const double targetFps = 60.0;
    const double framePeriod = 1.0 / targetFps;

    char title[128];

    while (!glfwWindowShouldClose(window)) {
        double now = glfwGetTime();

        // переключаем кадр раз в framePeriod секунд
        if (now - lastSwitch >= framePeriod) {
            currentFrame++;
            if (currentFrame >= numFrames) {
                currentFrame = numFrames - 1; // тормозим на последнем
                // если хочешь цикл по кругу:
                // currentFrame = 0;
            }

            char dsetName[64];
            std::snprintf(dsetName, sizeof(dsetName),
                          "frame_%06zu", currentFrame);

            try {
                auto parts = Reader(fileName, dsetName);
                if (parts.size() != objs.size()) {
                    std::cerr << "Размер кадра " << dsetName
                              << " не совпадает с количеством объектов\n";
                    break;
                }

                for (std::size_t i = 0; i < objs.size(); ++i) {
                    objs[i].position = parts[i].position;
                    objs[i].velocity = parts[i].velocity;
                    if (std::isfinite(parts[i].radius) && parts[i].radius > 0.0) {
                        objs[i].radius = static_cast<float>(parts[i].radius);
                    }
                    objs[i].UpdateVertices();
                }
            }
            catch (const std::exception& e) {
                std::cerr << "Ошибка при чтении " << dsetName << ": "
                          << e.what() << "\n";
                break;
            }

            lastSwitch = now;
        }

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        renderer.updateView(cameraPos, cameraFront, cameraUp);
        renderer.drawObjects(objs);

        std::snprintf(title, sizeof(title),
                      "Replay | frame %zu / %zu",
                      currentFrame + 1, numFrames);
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
