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
#include <filesystem>
#include "object.hpp"
#include "renderer.hpp"
#include "control.hpp"
#include "bodysystem.hpp"
#include "data.hpp"
#include "barnes_hut.hpp"
#include "H5Cpp.h"
#include <time.h>
#include "brut_force.hpp"
#include "physics.hpp"
#include "generators.hpp"
#include "camera.hpp"
#include "state.hpp"

#ifdef USE_CUDA
#include "barnes_hut_cuda.cuh"
#endif

#ifndef VEL_SCALE
#define VEL_SCALE 1.0f
#endif

#ifdef _WIN32
#include <windows.h>
#endif

float lastX = 400.0f, lastY = 300.0f;
float yaw   = -90.0f;
float pitch = 0.0f;

double gSimTime = 0.0;
float fixedDt = 1.0f/60; // шаг времени;
const int FIXED_STEPS = 10;
float initMass = 5.0f * std::pow(10.0f, 20.0f) / 5.0f;
char title[128];

std::vector<Object> objs = {};


GLFWwindow* StartGLU();


int main() {

#ifdef _WIN32
    SetConsoleOutputCP(65001);
    SetConsoleCP(65001);
#endif
    bool fullscreen = false;
    bool maximized = true;
    GLFWwindow* window = InitWindow(1280, 720, "3D_TEST", fullscreen, maximized);
    if (!window) {
        std::cerr << "Window or OpenGL context creation failed.\n";
        return -1;
    }

    int width, height;
    glfwGetFramebufferSize(window, &width, &height);

    Renderer renderer(width, height);
    renderer.setProjection(65.0f, (float) width/ (float) height, 8.3f, 100000.0f);
    using Handler = void(*)(std::vector<Object>& objs, float dt, bool pause, int iterations);
    Handler simulationStep = nullptr;

    RenderMode renderMode = RenderMode::Sphere;
    auto scenarioManager = CreateDefaultManager();

    bool loaded = false;
    int mode;
    std::cout << "Кубы/шары/точка [0/1/2]\n";
    std::cin >> mode;
    switch (mode){
        case 0:
            renderMode = RenderMode::Cubes;
            break;
        case 1:
            renderMode = RenderMode::Sphere;
            break;
        case 2:
            renderMode = RenderMode::Points;
            break;
        default:
            renderMode = RenderMode::Sphere;
            break;
    }
    renderer.setRenderMode(renderMode);

#ifdef USE_CUDA
    std::cout << "Выберите алгоритм:\n"
              << "  [1] Брутфорс (O(N^2))\n"
              << "  [2] Барнс-Хат CPU\n"
              << "  [3] Барнс-Хат GPU (CUDA)\n"
              << "Введите номер: " << std::flush;
#else
    std::cout << "Выберите алгоритм:\n"
              << "  [1] Брутфорс (O(N^2))\n"
              << "  [2] Барнс-Хат CPU\n"
              << "Введите номер: " << std::flush;
#endif
    std::cin >> mode;

    switch(mode) {
        case 1:
            simulationStep = &simulationStepBrutForceCPU;
            break;
        case 2:
            simulationStep = &simulationStepBarnesHutCPU;
            break;
#ifdef USE_CUDA
        case 3:
            simulationStep = &simulationStepBarnesHutCUDA;
            break;
#endif
        default:
            std::cout << "Неверный выбор, используем брутфорс.\n";
            simulationStep = &simulationStepBrutForceCPU;
            break;
    }

    // -------- выбор сценария (файл / рандом) --------
    std::cout << "Загружаем сценарий из HDF5 или генерируем систему рандомно? [0/1]: " << std::flush;
    std::cin >> mode;

    if (mode == 0){
        auto files = ListH5Files("data");
        if (files.empty()){
            std::cout << "В ./data нет .h5 файлов. Генерируем рандомно.\n";
        }
        else{
            for (size_t i = 0; i < files.size(); ++i){
                std::cout << "  [" << (i+1) << "] " << files[i] << "\n";
            }
            std::cout << "Выберите номер файла: " << std::flush;
            size_t idx = 0;
            std::cin >> idx;
            if (LoadObjectsFromFile(files[idx-1], "Particles", objs)) {
                physics::colorFromMass(objs);
                loaded = true;
                std::cout << "Loaded " << objs.size() << " objects\n";
            }

        }
    }
if (!loaded) {
        GenParams params;
        params.seed = 42;
        
        std::cout << "Количество тел: ";
        std::cin >> params.count;

        params.centralMass = static_cast<double>(initMass) * 100000;
        params.baseMass    = static_cast<double>(initMass);
        params.minRadius   = 300.0f;
        params.maxRadius   = 10000.0f;
        params.spread      = 1000.0f; 

        std::vector<std::string> scenNames = scenarioManager->getNames();
        std::cout << "Выберите тип генерации:\n";
        for (size_t i = 0; i < scenNames.size(); ++i) {
            std::cout << "  [" << (i + 1) << "] " << scenNames[i] << "\n";
        }
        
        size_t choice = 0;
        std::cout << "Ваш выбор: ";
        std::cin >> choice;
        
        if (choice > 0 && choice <= scenNames.size()) {
            scenarioManager->runScenario(choice - 1, objs, params);
        } else {
            std::cout << "Неверный выбор, запускаем первый по списку.\n";
            scenarioManager->runScenario(0, objs, params);
        }
    }

    H5::H5File framesFile = CreateSimulationFile("data/frames.h5", objs.size(), fixedDt * static_cast<double>(FIXED_STEPS));
    std::size_t frameIndex = 0;
    WriteSimulationFrame(framesFile, objs, frameIndex);
    ++frameIndex;
    // Чтение с HDF5

    BodySystem bodySystem(objs); // создание сохрянем информацию о системе

    bodySystem.transPointToSystem(objs); // переходим в систему объектов

    Camera cam;
    SimState state; 

    // Управление
    Control control(window, objs,
                    cam, state);
    control.attach();
    int counter = 0;
    double lastTime = glfwGetTime();
    double accumulator = 0.0;

    bool isRealTime;
    std::cout << "Режим реального времени? [1 - Да / 0 Нет]: " << std::flush;
    std::cin >> isRealTime;

    int stepCounter = 0;

    if (isRealTime) {
        while (!glfwWindowShouldClose(window) && state.running) {
            double now = glfwGetTime();
            double frameRealDt = now - lastTime;
            

            state.deltaTime = frameRealDt;
            lastTime = now;
            frameRealDt *= state.timeScale;
            accumulator += frameRealDt;

            int substeps = 0;

            while (accumulator >= fixedDt) {
                simulationStep(objs, fixedDt, state.pause, 1);
                gSimTime += fixedDt;
                accumulator -= fixedDt;
                ++substeps;
            }

            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            renderer.updateView(cam);
            renderer.drawObjects(objs);

            double dtForFps = frameRealDt / std::max(1.0f, state.timeScale);
            double fps = (dtForFps > 0.0) ? 1.0 / dtForFps : 0.0;
            std::snprintf(title, sizeof(title),
                  "REAL-TIME | Speed: %.1fx | FPS: %.0f | Obj: %zu | Time: %.2f",
                  state.timeScale, fps, objs.size(), gSimTime);
            glfwSetWindowTitle(window, title);

            glfwSwapBuffers(window);
            glfwPollEvents();
        }

    } else {
        int startTime = static_cast<int>(time(NULL));
        double targetTime;
        std::cout << "На какое время просчитываем? (сек): " << std::flush;
        std::cin >> targetTime;

        glfwSwapInterval(0);
        std::cout << "Начинаем расчет..." << std::endl;

        while (!glfwWindowShouldClose(window) && state.running && gSimTime < targetTime) {
            simulationStep(objs, fixedDt, false, FIXED_STEPS);
            gSimTime += fixedDt * FIXED_STEPS;
            stepCounter += FIXED_STEPS;
            WriteSimulationFrame(framesFile, objs, frameIndex);
            ++frameIndex;
            if (stepCounter % (FIXED_STEPS * 10) == 0) {
                glfwPollEvents();

                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
                renderer.updateView(cam);
                
                double progress = (gSimTime / targetTime) * 100.0;
                std::snprintf(title, sizeof(title),
                      "BAKING... %.1f%% | Time: %.2f / %.2f | Saved: %zu",
                      progress, gSimTime, targetTime, frameIndex);
                glfwSetWindowTitle(window, title);
                glfwSwapBuffers(window);
                
            }
        }

        std::cout << "Расчет завершен за " << (static_cast<int>(time(NULL)) - startTime) << " сек.\n";
    }
    CloseSimulationFile(framesFile, frameIndex);
}

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
