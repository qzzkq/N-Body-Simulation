#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <iostream>
#include <optional>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <filesystem>
#include <csignal>
#include <array>
#include <string>
#include <cstring>
// Добавлено для замера времени без GLFW
#include <chrono>
// Добавлено для std::unique_ptr
#include <memory>
#include "object.hpp"
#include "graphic_state.hpp"
#include "renderer.hpp"
#include "control.hpp"
#include "bodysystem.hpp"
#include "data.hpp"
#include "barnes_hut.hpp"
#include "H5Cpp.h"
#include <time.h>
#include "brut_force.hpp"
#include "whfast.hpp"
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

// Функции для парсинга аргументов командной строки 
std::string getCmdOption(int argc, char* argv[], const std::string& option, const std::string& default_val = "") {
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == option && i + 1 < argc) {
            return argv[i + 1];
        }
    }
    return default_val;
}

bool cmdOptionExists(int argc, char* argv[], const std::string& option) {
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == option) return true;
    }
    return false;
}

double gSimTime = 0.0;
// шаг времени
double fixedDt = 1.0 / 3652.5;
float initMass = 5.0f * std::pow(10.0f, 20.0f) / 5.0f;
char title[128];

std::vector<Object> objs = {};
std::vector<GraphicState> graphics;
volatile sig_atomic_t g_Interrupt = 0;

void signalHandler(int signum) {
    if (signum == SIGINT) {
        g_Interrupt = 1;
    }
}

int main(int argc, char* argv[]) {

    std::signal(SIGINT, signalHandler);

#ifdef _WIN32
    SetConsoleOutputCP(65001);
    SetConsoleCP(65001);
#endif

    bool is_cli = (argc > 1);

    bool isRealTime;
    if (!is_cli) {
        std::cout << "Режим реального времени (визуализация)? [1 - Да / 0 - Нет]: " << std::flush;
        std::cin >> isRealTime;
    } else {
        isRealTime = (getCmdOption(argc, argv, "--realtime", "1") == "1");
    }

    // --- Качество (скорость vs точность) ---
    // Сначала пресет (если задан), потом индивидуальные оверрайды.
    if (!is_cli) {
        std::cout << "Качество расчёта:\n"
                  << "  [1] Быстро (viz, ×3-5 быстрее, заметный энергодрейф)\n"
                  << "  [2] Сбалансированно (по умолчанию)\n"
                  << "  [3] Точно (медленно, минимальная погрешность)\n"
                  << "Ваш выбор: " << std::flush;
        int q = 2;
        std::cin >> q;
        physics::applyQualityPreset(std::clamp(q, 1, 3) - 1);
    } else {
        const std::string qStr = getCmdOption(argc, argv, "--quality", "");
        if (!qStr.empty()) {
            int mode = 1; // balanced
            if      (qStr == "fast"     || qStr == "1") mode = 0;
            else if (qStr == "balanced" || qStr == "2") mode = 1;
            else if (qStr == "precise"  || qStr == "3") mode = 2;
            else {
                std::cerr << "Неверное значение --quality (допустимо fast/balanced/precise), используем balanced\n";
            }
            physics::applyQualityPreset(mode);
        }
        const std::string thetaStr = getCmdOption(argc, argv, "--bh-theta", "");
        if (!thetaStr.empty()) {
            try {
                physics::setBarnesHutTheta(std::stod(thetaStr));
            } catch (...) {
                std::cerr << "Неверное значение --bh-theta, используем по умолчанию\n";
            }
        }
        const std::string softStr = getCmdOption(argc, argv, "--softening", "");
        if (!softStr.empty()) {
            try {
                physics::setSofteningAU(std::stod(softStr));
            } catch (...) {
                std::cerr << "Неверное значение --softening, используем по умолчанию\n";
            }
        }
        const std::string etaStr = getCmdOption(argc, argv, "--adaptive-eta", "");
        if (!etaStr.empty()) {
            try {
                physics::setAdaptiveEta(std::stod(etaStr));
            } catch (...) {
                std::cerr << "Неверное значение --adaptive-eta, используем по умолчанию\n";
            }
        }
        const std::string capStr = getCmdOption(argc, argv, "--substep-cap", "");
        if (!capStr.empty()) {
            try {
                physics::setSubstepCap(std::stoi(capStr));
            } catch (...) {
                std::cerr << "Неверное значение --substep-cap, используем по умолчанию\n";
            }
        }
    }
    std::cout << "Параметры: theta=" << physics::getBarnesHutTheta()
              << ", softening=" << physics::getSofteningAU()
              << ", eta=" << physics::getAdaptiveEta()
              << ", cap=" << physics::getSubstepCap() << "\n";

    Renderer renderer; 
    std::unique_ptr<Control> control = nullptr;
    RenderMode renderMode = RenderMode::Sphere;

    if (isRealTime) {
        bool fullscreen = false;
        bool maximized = true;

        if (!renderer.init(1280, 720, "N-Body simulation", fullscreen, maximized)) {
            std::cerr << "Window or OpenGL context creation failed.\n";
            return -1;
        }

        renderer.setProjection(65.0f, 1.0f, 1.0e10f);

        int mode;
        if (!is_cli) {
            std::cout << "Кубы/шары/точка [0/1/2]\n";
            std::cin >> mode;
        } else {
            mode = std::stoi(getCmdOption(argc, argv, "--render", "1"));
        }

        switch (mode){
            case 0: renderMode = RenderMode::Cubes; break;
            case 1: renderMode = RenderMode::Sphere; break;
            case 2: renderMode = RenderMode::Points; break;
            default: renderMode = RenderMode::Sphere; break;
        }
        renderer.setRenderMode(renderMode);
    }

    using Handler = void(*)(std::vector<Object>& objs, double dt, bool pause, bool forceSync);
    Handler simulationStep = nullptr;
    auto scenarioManager = CreateDefaultManager();

    bool loaded = false;
    bool loadedColorsFromTxt = false;
    int mode;

    // --- 2. Выбор алгоритма ---
#ifdef USE_CUDA
    if (!is_cli) {
        std::cout << "Выберите алгоритм:\n"
                  << "  [1] Брутфорс Yoshida4 (O(N^2))\n"
                  << "  [2] Барнс-Хат CPU\n"
                  << "  [3] Барнс-Хат GPU (CUDA)\n"
                  << "  [4] WHFast (Solar-like, O(N^2) на kick + Kepler)\n"
                  << "Введите номер: " << std::flush;
        std::cin >> mode;
    } else {
        mode = std::stoi(getCmdOption(argc, argv, "--algo", "1"));
    }
#else
    if (!is_cli) {
        std::cout << "Выберите алгоритм:\n"
                  << "  [1] Брутфорс Yoshida4 (O(N^2))\n"
                  << "  [2] Барнс-Хат CPU\n"
                  << "  [4] WHFast (Solar-like, O(N^2) на kick + Kepler)\n"
                  << "Введите номер: " << std::flush;
        std::cin >> mode;
    } else {
        mode = std::stoi(getCmdOption(argc, argv, "--algo", "1"));
    }
#endif

    switch(mode) {
        case 1: simulationStep = &simulationStepBrutForceCPU; break;
        case 2: simulationStep = &simulationStepBarnesHutCPU; break;
#ifdef USE_CUDA
        case 3: simulationStep = &simulationStepBarnesHutCUDA; break;
#endif
        case 4: simulationStep = &simulationStepWHFastCPU; break;
        default:
            std::cout << "Неверный выбор, используем брутфорс.\n";
            simulationStep = &simulationStepBrutForceCPU;
            break;
    }

    // --- 3. Шаг времени dt ---
    std::string dtInput;
    if (!is_cli) {
        std::cout << "Введите шаг времени dt (например, 0.000273785 для 1/3652.5): ";
        std::cin >> dtInput;
    } else {
        dtInput = getCmdOption(argc, argv, "--dt", "0.000273785");
    }
    
    try {
        fixedDt = std::stod(dtInput);
    } catch (...) {
        std::cout << "Некорректный ввод. Используется dt по умолчанию (1/3652.5).\n";
        fixedDt = 1.0 / 3652.5;
    }

    // --- 4. Источник начальных условий ---
    if (!is_cli) {
        std::cout << "Источник начальных условий: [0] HDF5, [1] TXT, [2] Random: " << std::flush;
        std::cin >> mode;
    } else {
        mode = std::stoi(getCmdOption(argc, argv, "--source", "2"));
    }

    if (mode == 0) {
        if (!is_cli) {
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
                bool h5HasColors = false;
                if (LoadObjectsFromFile(files[idx-1], "Particles", objs, &graphics, &h5HasColors)) {
                    loaded = true;
                    if (h5HasColors) loadedColorsFromTxt = true;
                    std::cout << "Loaded " << objs.size() << " objects\n";
                }
            }
        } else {
            std::string h5path = getCmdOption(argc, argv, "--h5", "");
            bool h5HasColors = false;
            if (!h5path.empty() && LoadObjectsFromFile(h5path, "Particles", objs, &graphics, &h5HasColors)) {
                loaded = true;
                if (h5HasColors) loadedColorsFromTxt = true;
                std::cout << "Loaded " << objs.size() << " objects\n";
            } else {
                std::cout << "Failed to load HDF5 from path: " << h5path << "\n";
            }
        }
    }
    else if (mode == 1) {
        std::string txtPath;
        if (!is_cli) {
            std::cout << "Введите путь к TXT-файлу (пример: data/solar_system.txt): " << std::flush;
            std::cin >> txtPath;
        } else {
            txtPath = getCmdOption(argc, argv, "--txt", "");
        }

        if (!txtPath.empty() && LoadSystemFromTextFile(txtPath, objs, &graphics)) {
            loaded = true;
            loadedColorsFromTxt = true;
            std::cout << "Loaded " << objs.size() << " objects from text config\n";
        } else {
            std::cout << "Не удалось загрузить TXT. Генерируем рандомно.\n";
        }
    }

    if (!loaded) {
        GenParams params;
        params.seed = 42;
        
        if (!is_cli) {
            std::cout << "Количество тел: ";
            std::cin >> params.count;
        } else {
            params.count = std::stoi(getCmdOption(argc, argv, "--bodies", "100"));
        }

        params.centralMass = 1.0;      
        params.baseMass    = 3.0e-6;   
        params.minRadius   = 5.0f;     
        params.maxRadius   = 40.0f;    
        params.spread      = 2.0f;     

        std::vector<std::string> scenNames = scenarioManager->getNames();
        size_t choice = 0;
        
        if (!is_cli) {
            std::cout << "Выберите тип генерации:\n";
            for (size_t i = 0; i < scenNames.size(); ++i) {
                std::cout << "  [" << (i + 1) << "] " << scenNames[i] << "\n";
            }
            std::cout << "Ваш выбор: ";
            std::cin >> choice;
        } else {
            choice = std::stoi(getCmdOption(argc, argv, "--scenario", "0")) + 1;
        }
        
        if (choice > 0 && choice <= scenNames.size()) {
            scenarioManager->runScenario(choice - 1, objs, params);
        } else {
            std::cout << "Неверный выбор, запускаем первый по списку.\n";
            scenarioManager->runScenario(0, objs, params);
        }
    }

    if (!loadedColorsFromTxt) {
        graphics.clear();
        graphics.resize(objs.size());
        physics::colorFromMass(objs, graphics);
    }

    // --- 5. Настройки сохранения ---
    std::string filename;
    if (!is_cli) {
        std::cout << "Укажите название файла для сохранения (без расширения): " << std::flush;
        std::cin >> filename;
    } else {
        filename = getCmdOption(argc, argv, "--output", "frames");
    }

    int saveIntervalSteps = 10;
    if (!is_cli) {
        std::cout << "Как часто сохранять кадры в HDF5? Целое N >= 1: каждый N-й шаг интеграции "
                     "(1 = каждый шаг, 10 = раз в 10 шагов).\nN = " << std::flush;
        std::cin >> saveIntervalSteps;
    } else {
        saveIntervalSteps = std::stoi(getCmdOption(argc, argv, "--save_every", "10"));
    }

    if (saveIntervalSteps < 1) {
        if (!is_cli) std::cout << "N < 1, используем N = 1.\n";
        saveIntervalSteps = 1;
    }

    std::filesystem::create_directories("data");
    H5::H5File framesFile = CreateSimulationFile("data/" + filename + ".h5", objs.size(),
                                               fixedDt * static_cast<double>(saveIntervalSteps));
    std::size_t frameIndex = 0;
    WriteSimulationFrame(framesFile, objs, graphics, frameIndex);
    ++frameIndex;
    
    BodySystem bodySystem(objs); 
    bodySystem.transPointToSystem(objs); 

    const double initialEnergy = physics::calculateTotalEnergy(objs);
    double energyCache = initialEnergy;

    Camera cam;
    SimState state; 

    int stepCounter = 0;
    int frameCounter = 0;

    // =========================================================
    // ЦИКЛ СИМУЛЯЦИИ (РЕАЛТАЙМ ИЛИ БЕЙКИНГ)
    // =========================================================
    if (isRealTime) {
        control = std::make_unique<Control>(renderer.getWindow(), objs, cam, state);
        control->attach();

        double lastTime = glfwGetTime();
        double accumulator = 0.0;

        while (!glfwWindowShouldClose(renderer.getWindow()) && state.running && !g_Interrupt) {
            double now = glfwGetTime();
            double frameRealDt = now - lastTime;
            
            state.deltaTime = frameRealDt;
            lastTime = now;
            frameRealDt *= state.timeScale;
            if (frameRealDt > 0.1) frameRealDt = 0.1;
            if (frameRealDt < 0.0) frameRealDt = 0.0;
            
            GLFWwindow* win = renderer.getWindow();
            const bool shiftFreeze = (glfwGetKey(win, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
                || (glfwGetKey(win, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS);
            if (!shiftFreeze) {
                accumulator += frameRealDt;
            }

            int substeps = 0;
            constexpr int MAX_SUBSTEPS = 100000;
            while (accumulator >= fixedDt && substeps < MAX_SUBSTEPS && !g_Interrupt) {
                // realtime: после каждого шага рендер читает objs[].position →
                // нужен host-sync.
                simulationStep(objs, fixedDt, state.pause, true);
                gSimTime += fixedDt;
                accumulator -= fixedDt;
                ++substeps;
                for (std::size_t i = 0; i < objs.size(); ++i) {
                    graphics[i].updateTrail(objs[i].position);
                }
            }
            if (substeps >= MAX_SUBSTEPS) {
                accumulator = 0.0;
            }

            control->updateCameraFromKeys();
            renderer.renderFrame(objs, graphics, cam);
            ++frameCounter;

            double dtForFps = frameRealDt / std::max(1.0f, state.timeScale);
            double fps = (dtForFps > 0.0) ? 1.0 / dtForFps : 0.0;
            
            std::snprintf(title, sizeof(title),
                  "REAL-TIME | Sim: %.1fx | Cam: %.2fx | FPS: %.0f | Obj: %zu | Time: %.2f%s",
                  state.timeScale, control->getCameraMoveScale(), fps, objs.size(), gSimTime,
                  shiftFreeze ? " | [SHIFT freeze]" : "");
            glfwSetWindowTitle(renderer.getWindow(), title);

            glfwPollEvents();
        }

    } else {
        int startTime = static_cast<int>(time(NULL));
        double targetTime = 0.0;
        
        {
            const double yearsBetweenSaves = fixedDt * static_cast<double>(saveIntervalSteps);
            const double mbPerSimYear = (yearsBetweenSaves > 0.0 && !objs.empty())
                    ? (static_cast<double>(objs.size()) * 3.0 * static_cast<double>(sizeof(float))
                       / (yearsBetweenSaves * 1024.0 * 1024.0)) : 0.0;
                    
            if (!is_cli) {
                std::cout << std::fixed << std::setprecision(6);
                std::cout << "\nПросчёт (не в реальном времени). Один шаг dt = " << fixedDt << " года.\n";
                std::cout << "Между сохранёнными кадрами: " << yearsBetweenSaves << " года ("
                          << saveIntervalSteps << " шаг(ов)). Ориентир размера: ~" << std::setprecision(3)
                          << mbPerSimYear << " МБ на 1 год симуляции.\n";
                std::cout << "На какое симуляционное время просчитываем? (годы): " << std::flush;
                std::cin >> targetTime;
            } else {
                targetTime = std::stod(getCmdOption(argc, argv, "--target", "10.0"));
            }
        }
        
        // Используем std::chrono вместо glfwGetTime()
        auto bakeStartRealTime = std::chrono::high_resolution_clock::now();
        std::cout << "Начинаем расчет (headless режим)..." << std::endl;

        while (state.running && gSimTime < targetTime && !g_Interrupt) {
            double stepDt = std::min(fixedDt, targetTime - gSimTime);
            // headless: host-sync только когда кадр будет сохранён.
            // Это убирает D→H копию (pos+vel+acc, ~5 МБ) на каждом не-save шаге.
            const bool willSave = (((stepCounter + 1) % saveIntervalSteps) == 0);
            simulationStep(objs, stepDt, false, willSave);
            gSimTime += stepDt;
            stepCounter += 1;
            
            if (stepCounter % (1 * 10) == 0) {
                double progress = (gSimTime / targetTime) * 100.0;
                auto now = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed_chrono = now - bakeStartRealTime;
                double elapsed = elapsed_chrono.count();
                
                double progress_dec = gSimTime / targetTime;
                double etaSeconds = 0.0;
                if (progress_dec > 0.001) {
                    etaSeconds = (elapsed / progress_dec) - elapsed;
                    if (etaSeconds < 0.0) etaSeconds = 0.0;
                }
                int etaTotalSec = static_cast<int>(etaSeconds);
                int etaMin = etaTotalSec / 60;
                int etaSec = etaTotalSec % 60;
                
                std::cout << "\rBAKING... " << std::fixed << std::setprecision(1) << progress 
                          << "% | ETA: " << std::setfill('0') << std::setw(2) << etaMin << ":" 
                          << std::setw(2) << etaSec 
                          << " | Time: " << std::setprecision(2) << gSimTime << " / " << targetTime 
                          << " | Saved: " << frameIndex << std::flush;
            }
            if (stepCounter % saveIntervalSteps == 0) {
                WriteSimulationFrame(framesFile, objs, graphics, frameIndex);
                ++frameIndex;
            }
        }
        std::cout << "\nРасчет завершен за " << (static_cast<int>(time(NULL)) - startTime) << " сек.\n";
    }

    if (g_Interrupt) {
        std::cout << "\nПолучен SIGINT (Ctrl+C). Завершаем симуляцию и сохраняем текущее состояние...\n";
    }

    simulationStep(objs, 0.0, true, true);
    const std::string finalTextPath = "data/" + filename + "_final.txt";
    if (SaveSystemToTextFile(finalTextPath, objs, &graphics)) {
        std::cout << "Финальное состояние сохранено в " << finalTextPath << "\n";
    } else {
        std::cout << "Не удалось сохранить финальное состояние в " << finalTextPath << "\n";
    }

    CloseSimulationFile(framesFile, frameIndex);
    return 0;
}