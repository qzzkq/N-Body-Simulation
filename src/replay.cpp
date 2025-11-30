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

// ===== тип первого кадра: position + mass + radius =====
struct Frame0Record {
    glm::dvec3 position;
    double     mass;
    double     radius;
};

static H5::CompType MakeFrame0TypeReplay()
{
    hsize_t v3dims[1] = {3};
    H5::ArrayType vec3Type(H5::PredType::NATIVE_DOUBLE, 1, v3dims);

    H5::CompType t(sizeof(Frame0Record));
    t.insertMember("position", HOFFSET(Frame0Record, position), vec3Type);
    t.insertMember("mass",     HOFFSET(Frame0Record, mass),     H5::PredType::NATIVE_DOUBLE);
    t.insertMember("radius",   HOFFSET(Frame0Record, radius),   H5::PredType::NATIVE_DOUBLE);
    return t;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Использование: " << argv[0] << " output.h5\n";
        return 1;
    }

    std::string fileName = argv[1];

    H5::H5File file(fileName, H5F_ACC_RDONLY);

    // === детектируем количество кадров по наличию frame_XXXXXX ===
    std::size_t numFrames = 0;
    while (true) {
        char dsetName[64];
        std::snprintf(dsetName, sizeof(dsetName), "frame_%06zu", numFrames);
        try {
            file.openDataSet(dsetName).close();
            ++numFrames;
        } catch (const H5::Exception&) {
            break;
        }
    }

    if (numFrames == 0) {
        std::cerr << "В файле " << fileName
                  << " не найдено ни одного кадра frame_XXXXXX.\n";
        return 1;
    }

    // === читаем первый кадр (position + mass + radius) ===
    char firstName[64];
    std::snprintf(firstName, sizeof(firstName), "frame_%06d", 0);

    H5::DataSet ds0 = file.openDataSet(firstName);
    H5::DataSpace sp0 = ds0.getSpace();

    hsize_t dims[1] = {0};
    sp0.getSimpleExtentDims(dims);
    std::size_t numBodies = static_cast<std::size_t>(dims[0]);
    if (numBodies == 0) {
        std::cerr << "Первый кадр пустой.\n";
        return 1;
    }

    std::vector<Frame0Record> firstFrame(numBodies);
    H5::CompType t0 = MakeFrame0TypeReplay();
    ds0.read(firstFrame.data(), t0);

    std::cout << "Файл: " << fileName
              << "\n  num_bodies (detected) = " << numBodies
              << "\n  num_frames (detected) = " << numFrames << "\n";

    // === OpenGL / окно ===
    GLFWwindow* window = StartGLU();
    if (!window) {
        return 1;
    }

    Renderer renderer(800, 600, vertexShaderSource, fragmentShaderSource);
    renderer.setProjection(65.0f, 800.0f / 600.0f, 8.3f, 100000.0f);

    glm::vec3 cameraPos   = glm::vec3(0.0f, 50.0f, 250.0f);
    glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
    glm::vec3 cameraUp    = glm::vec3(0.0f, 1.0f,  0.0f);

    // === создаём объекты по первому кадру ===
    std::vector<Object> objs;
    objs.reserve(numBodies);
    for (const auto& r : firstFrame) {
        Object o(r.position,
                 glm::dvec3(0.0),  // скорость в реплее не нужна
                 r.mass,
                 1410.0f,
                 std::nullopt);
        o.Initalizing = false;
        o.radius = static_cast<float>(r.radius);
        if (!std::isfinite(o.radius) || o.radius <= 0.0f) {
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

    // тип для последующих кадров (только position)
    hsize_t v3dims[1] = {3};
    H5::ArrayType vec3Type(H5::PredType::NATIVE_DOUBLE, 1, v3dims);

    while (!glfwWindowShouldClose(window)) {
        double now = glfwGetTime();

        if (now - lastSwitch >= framePeriod) {
            currentFrame++;
            if (currentFrame >= numFrames) {
                currentFrame = numFrames - 1; // стоп на последнем
                // или по кругу:
                // currentFrame = 0;
            }

            char dsetName[64];
            std::snprintf(dsetName, sizeof(dsetName),
                          "frame_%06zu", currentFrame);

            try {
                H5::DataSet dset = file.openDataSet(dsetName);
                H5::DataSpace sp = dset.getSpace();

                hsize_t dimsPos[1] = {0};
                sp.getSimpleExtentDims(dimsPos);
                if (dimsPos[0] != numBodies) {
                    std::cerr << "Размер кадра " << dsetName
                              << " не совпадает с numBodies\n";
                    break;
                }

                std::vector<glm::dvec3> pos(numBodies);
                dset.read(pos.data(), vec3Type);

                for (std::size_t i = 0; i < objs.size(); ++i) {
                    objs[i].position = pos[i];
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
