#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <H5Cpp.h>

#include <iostream>
#include <vector>
#include <string>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <limits>

#ifdef _WIN32
#include <windows.h>
#endif

#include "object.hpp"
#include "renderer.hpp"
#include "control.hpp"

bool running = true;
bool pause   = false;

glm::vec3 cameraPos   = glm::vec3(0.0f, 0.0f, 1.0f);
glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
glm::vec3 cameraUp    = glm::vec3(0.0f, 1.0f,  0.0f);

float lastX = 400.0f, lastY = 300.0f;
float yaw   = -90.0f;
float pitch = 0.0f;

float dt        = 0.0f;
float timeScale = 1.0f;   // скорость воспроизведения

char title[256];

std::vector<Object> objs;


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

// Инициализируем контекст openGL
GLFWwindow* StartGLU() {

    if (!glfwInit()) {
        std::cout << "Failed to initialize GLFW\n";
        return nullptr;
    }

    // получаем параметры монитора 
    GLFWmonitor* monitor = glfwGetPrimaryMonitor(); 
    if (monitor == NULL) {
        std::cerr << "Failed to create GLFW window";
        return nullptr; 
    }

    const GLFWvidmode* mode = glfwGetVideoMode(monitor);
    if (mode == NULL) {
        std::cerr << "Failed to create GLFW window";
        return nullptr; 
    }

    GLFWwindow* window = glfwCreateWindow(mode->width, mode->height, "3D_TEST", nullptr, nullptr);
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
    glViewport(0, 0, mode->width, mode->height);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    return window;
}


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

void colorFromMass(std::vector<Object>& objs) {
    if (objs.empty()) return;

    double minMass = std::numeric_limits<double>::max();
    double maxMass = -1.0;

    // определение минимальной массы 
    for (const auto& obj : objs) {
        if (obj.mass > maxMass) maxMass = obj.mass;
        if (obj.mass < minMass) minMass = obj.mass;
    }

    // Защита от нулевой массы и log(0)
    minMass = std::max(minMass, 1e-30);
    maxMass = std::max(maxMass, 1e-30);

    // 2. Палитра цветов для планет 
    const glm::vec3 colorCold = glm::vec3(0.8f, 0.4f, 0.0f); // Темно-оранжевый (лёгкие)
    const glm::vec3 colorMid  = glm::vec3(1.0f, 0.9f, 0.1f); // Ярко-желтый (средние и тяжелее среднего)
    const glm::vec3 colorHot  = glm::vec3(1.0f, 1.0f, 1.0f); // Белый (самые тяжёлые)

    // Используем логарифмическую шкалу, так как массы планет отличаются в миллиарды раз
    double logMin = std::log10(minMass);
    double logMax = std::log10(maxMass);
    double range = logMax - logMin;

    // присваиваем цвета объектам
    for (auto& obj : objs) {
        double m = std::max(obj.mass, 1e-30);
        double val = std::log10(m);
        float t = 0.0f;
        if (range > 1e-9) {
            t = static_cast<float>((val - logMin) / range);
        }
        t = std::clamp(t, 0.0f, 1.0f);
        glm::vec3 finalColor;
        if (t < 0.6f) {
            float localT = t / 0.6f; 
            finalColor = glm::mix(colorCold, colorMid, localT);
        } 
        else {
            float localT = (t - 0.6f) / 0.4f;
            finalColor = glm::mix(colorMid, colorHot, localT);
        }
        obj.color = glm::vec4(finalColor, 1.0f);
    }
}

int main(int argc, char** argv) {
#ifdef _WIN32
    SetConsoleOutputCP(65001);
    SetConsoleCP(65001);
#endif

    H5::Exception::dontPrint();

    if (argc < 2) {
        std::cerr << "Использование: " << argv[0] << " data/frames.h5\n";
        return 1;
    }

    std::string fileName = argv[1];

    H5::H5File file(fileName, H5F_ACC_RDONLY);

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

    double initialTime = 0.0;
    try {
        H5::Attribute timeAttr = ds0.openAttribute("time");
        timeAttr.read(H5::PredType::NATIVE_DOUBLE, &initialTime);
    } catch (...) {
        initialTime = 0.0;
    }

    std::cout << "Файл: " << fileName
              << "\n  num_bodies (detected) = " << numBodies
              << "\n  num_frames (detected) = " << numFrames
              << "\n  t0 = " << initialTime << " s\n";

    GLFWwindow* window = StartGLU();
    if (!window) {
        return 1;
    }

    int width, height;
    glfwGetFramebufferSize(window, &width, &height);

    Renderer renderer(width, height, vertexShaderSource, fragmentShaderSource);
    renderer.setProjection(65.0f, (float) width/ (float) height, 8.3f, 100000.0f);

    cameraPos   = glm::vec3(0.0f, 50.0f, 250.0f);
    cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
    cameraUp    = glm::vec3(0.0f, 1.0f,  0.0f);

    objs.clear();
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

    colorFromMass(objs);

    float initMass = 5.0f * std::pow(10.0f, 20.0f) / 5.0f;

    Control control(window, objs,
                    cameraPos, cameraFront, cameraUp,
                    dt, timeScale, pause, running,
                    yaw, pitch, lastX, lastY,
                    initMass);
    control.attach();

    std::size_t currentFrame = 0;
    double currentSimTime    = initialTime;

    hsize_t v3dims[1] = {3};
    H5::ArrayType vec3Type(H5::PredType::NATIVE_DOUBLE, 1, v3dims);

    const double baseFrameDt = 1.0 / 60.0; // виртуальный шаг между кадрами
    double frameAccumulator  = 0.0;

    double lastTime = glfwGetTime();

    while (!glfwWindowShouldClose(window) && running) {
        double now = glfwGetTime();
        double realDt = now - lastTime;
        lastTime = now;

        dt = static_cast<float>(realDt);

        float clampedTimeScale = std::max(0.0f, timeScale);

        if (!pause && clampedTimeScale > 0.0f && numFrames > 1) {
            frameAccumulator += realDt * static_cast<double>(clampedTimeScale);

            while (frameAccumulator >= baseFrameDt &&
                   currentFrame + 1 < numFrames) {
                ++currentFrame;
                frameAccumulator -= baseFrameDt;

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

                    currentSimTime = currentFrame * baseFrameDt;
                    try {
                        H5::Attribute timeAttr = dset.openAttribute("time");
                        timeAttr.read(H5::PredType::NATIVE_DOUBLE, &currentSimTime);
                    } catch (...) {
                    }
                }
                catch (const std::exception& e) {
                    std::cerr << "Ошибка при чтении " << dsetName << ": "
                              << e.what() << "\n";
                    break;
                }
            }
        }

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        renderer.updateView(cameraPos, cameraFront, cameraUp);
        renderer.drawObjects(objs);

        std::snprintf(title, sizeof(title),
                      "Replay | frame %zu / %zu | t = %.3f s | speed x%.2f",
                      currentFrame, numFrames - 1,
                      currentSimTime, clampedTimeScale);
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
