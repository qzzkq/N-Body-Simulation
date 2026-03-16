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
#include <csignal>
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

double gSimTime = 0.0;
double fixedDt = 1.0 / 3652.5; // шаг времени;
const int FIXED_STEPS = 10;
float initMass = 5.0f * std::pow(10.0f, 20.0f) / 5.0f;
char title[128];

std::vector<Object> objs = {};
volatile sig_atomic_t g_Interrupt = 0;

void signalHandler(int signum) {
    if (signum == SIGINT) {
        g_Interrupt = 1;
    }
}

int main() {

    std::signal(SIGINT, signalHandler);

#ifdef _WIN32
    SetConsoleOutputCP(65001);
    SetConsoleCP(65001);
#endif
    bool fullscreen = false;
    bool maximized = true;

    Renderer renderer; 

    if (!renderer.init(1280, 720, "N-Body simulation", fullscreen, maximized)) {
        std::cerr << "Window or OpenGL context creation failed.\n";
        return -1;
    }

    renderer.setProjection(65.0f, 0.01f, 500.0f);
    using Handler = void(*)(std::vector<Object>& objs, double dt, bool pause);
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

    std::cout << "Введите шаг времени dt (например, 0.000273785 для 1/3652.5): ";
    std::string dtInput;
    std::cin >> dtInput;
    
    try {
        fixedDt = std::stod(dtInput);
    } catch (...) {
        std::cout << "Некорректный ввод. Используется dt по умолчанию (1/3652.5).\n";
        fixedDt = 1.0 / 3652.5;
    }

    // -------- выбор сценария (HDF5 / TXT / рандом) --------
    std::cout << "Источник начальных условий: [0] HDF5, [1] TXT, [2] Random: " << std::flush;
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
    else if (mode == 1) {
        std::string txtPath;
        std::cout << "Введите путь к TXT-файлу (пример: data/solar_system.txt): " << std::flush;
        std::cin >> txtPath;

        if (LoadSystemFromTextFile(txtPath, objs)) {
            loaded = true;
            std::cout << "Loaded " << objs.size() << " objects from text config\n";
        } else {
            std::cout << "Не удалось загрузить TXT. Генерируем рандомно.\n";
        }
    }

if (!loaded) {
        GenParams params;
        params.seed = 42;
        
        std::cout << "Количество тел: ";
        std::cin >> params.count;

        params.centralMass = 1.0;      // 1 масса Солнца
        params.baseMass    = 3.0e-6;   // Масса Земли в солнечных
        params.minRadius   = 5.0f;     // AU
        params.maxRadius   = 40.0f;    // AU
        params.spread      = 2.0f;     // AU

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

    std::cout << "Укажите название файла для сохранения (без расширения): " << std::flush;
    std::string filename;
    std::cin >> filename;
    H5::H5File framesFile = CreateSimulationFile("data/" + filename + ".h5", objs.size(), fixedDt * FIXED_STEPS);
    std::size_t frameIndex = 0;
    WriteSimulationFrame(framesFile, objs, frameIndex);
    ++frameIndex;
    // Чтение с HDF5

    BodySystem bodySystem(objs); // создание сохрянем информацию о системе

    bodySystem.transPointToSystem(objs); // переходим в систему объектов

    const double initialEnergy = physics::calculateTotalEnergy(objs);

    Camera cam;
    SimState state; 

    // Управление
    Control control(renderer.getWindow(), objs,
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
        while (!glfwWindowShouldClose(renderer.getWindow()) && state.running && !g_Interrupt) {
            double now = glfwGetTime();
            double frameRealDt = now - lastTime;
            

            state.deltaTime = frameRealDt;
            lastTime = now;
            frameRealDt *= state.timeScale;
            accumulator += frameRealDt;

            int substeps = 0;

            while (accumulator >= fixedDt && !g_Interrupt) {
                simulationStep(objs, fixedDt, state.pause);
                gSimTime += fixedDt;
                accumulator -= fixedDt;
                ++substeps;
                for(auto& obj : objs) obj.updateTrail();
            }

            //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            //renderer.updateView(cam);
            //renderer.drawObjects(objs);

            renderer.renderFrame(objs, cam);
            double dtForFps = frameRealDt / std::max(1.0f, state.timeScale);
            double fps = (dtForFps > 0.0) ? 1.0 / dtForFps : 0.0;
            const double currentEnergy = physics::calculateTotalEnergy(objs);
            const double energyDenom = std::max(std::abs(initialEnergy), 1e-30);
            const double errorPct = std::abs((currentEnergy - initialEnergy) / energyDenom) * 100.0;
            std::snprintf(title, sizeof(title),
                  "REAL-TIME | Speed: %.1fx | FPS: %.0f | Obj: %zu | Time: %.2f | Err: %.4f%%",
                  state.timeScale, fps, objs.size(), gSimTime, errorPct);
            glfwSetWindowTitle(renderer.getWindow(), title);

            //glfwSwapBuffers(renderer.getWindow());
            glfwPollEvents();
        }

    } else {
        int startTime = static_cast<int>(time(NULL));
        double targetTime;
        std::cout << "На какое время просчитываем? (сек): " << std::flush;
        std::cin >> targetTime;

        glfwSwapInterval(0);
        std::cout << "Начинаем расчет..." << std::endl;

        while (!glfwWindowShouldClose(renderer.getWindow()) && state.running && gSimTime < targetTime && !g_Interrupt) {
            simulationStep(objs, fixedDt, false);
            gSimTime += fixedDt;
            stepCounter += 1;
            if (stepCounter % (1 * 10) == 0) {
                glfwPollEvents();

                renderer.renderFrame(objs, cam);
                double progress = (gSimTime / targetTime) * 100.0;
                const double currentEnergy = physics::calculateTotalEnergy(objs);
                const double energyDenom = std::max(std::abs(initialEnergy), 1e-30);
                const double errorPct = std::abs((currentEnergy - initialEnergy) / energyDenom) * 100.0;
                std::snprintf(title, sizeof(title),
                      "BAKING... %.1f%% | Time: %.2f / %.2f | Saved: %zu | Err: %.4f%%",
                      progress, gSimTime, targetTime, frameIndex, errorPct);
                glfwSetWindowTitle(renderer.getWindow(), title);
            }
            if (stepCounter % FIXED_STEPS == 0) {
                WriteSimulationFrame(framesFile, objs, frameIndex);
                ++frameIndex;
            }
        }
        std::cout << "Расчет завершен за " << (static_cast<int>(time(NULL)) - startTime) << " сек.\n";
    }

    if (g_Interrupt) {
        std::cout << "Получен SIGINT (Ctrl+C). Завершаем симуляцию и сохраняем текущее состояние...\n";
    }

    const std::string finalTextPath = "data/" + filename + "_final.txt";
    if (SaveSystemToTextFile(finalTextPath, objs)) {
        std::cout << "Финальное состояние сохранено в " << finalTextPath << "\n";
    } else {
        std::cout << "Не удалось сохранить финальное состояние в " << finalTextPath << "\n";
    }

    CloseSimulationFile(framesFile, frameIndex);
}
