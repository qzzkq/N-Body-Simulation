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
double fixedDt = 1.0 / 3652.5; // шаг времени;
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

    bool is_cli = (argc > 1); // Флаг: запущен ли с аргументами командной строки

    bool fullscreen = false;
    bool maximized = true;

    Renderer renderer; 

    if (!renderer.init(1280, 720, "N-Body simulation", fullscreen, maximized)) {
        std::cerr << "Window or OpenGL context creation failed.\n";
        return -1;
    }

    // Дальняя плоскость в AU (логарифмический Z в шейдерах): ~10^10 AU — очень большая дальность видения.
    renderer.setProjection(65.0f, 1.0f, 1.0e10f);
    using Handler = void(*)(std::vector<Object>& objs, double dt, bool pause, bool forceSync);
    Handler simulationStep = nullptr;

    RenderMode renderMode = RenderMode::Sphere;
    auto scenarioManager = CreateDefaultManager();

    bool loaded = false;
    bool loadedColorsFromTxt = false;
    int mode;
    
    // --- 1. Выбор рендера ---
    if (!is_cli) {
        std::cout << "Кубы/шары/точка [0/1/2]\n";
        std::cin >> mode;
    } else {
        mode = std::stoi(getCmdOption(argc, argv, "--render", "1"));
    }

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

    // --- 2. Выбор алгоритма ---
#ifdef USE_CUDA
    if (!is_cli) {
        std::cout << "Выберите алгоритм:\n"
                  << "  [1] Брутфорс (O(N^2))\n"
                  << "  [2] Барнс-Хат CPU\n"
                  << "  [3] Барнс-Хат GPU (CUDA)\n"
                  << "Введите номер: " << std::flush;
        std::cin >> mode;
    } else {
        mode = std::stoi(getCmdOption(argc, argv, "--algo", "1"));
    }
#else
    if (!is_cli) {
        std::cout << "Выберите алгоритм:\n"
                  << "  [1] Брутфорс (O(N^2))\n"
                  << "  [2] Барнс-Хат CPU\n"
                  << "Введите номер: " << std::flush;
        std::cin >> mode;
    } else {
        mode = std::stoi(getCmdOption(argc, argv, "--algo", "1"));
    }
#endif

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
                    if (h5HasColors) {
                        loadedColorsFromTxt = true;
                    }
                    std::cout << "Loaded " << objs.size() << " objects\n";
                }
            }
        } else {
            std::string h5path = getCmdOption(argc, argv, "--h5", "");
            bool h5HasColors = false;
            if (!h5path.empty() && LoadObjectsFromFile(h5path, "Particles", objs, &graphics, &h5HasColors)) {
                loaded = true;
                if (h5HasColors) {
                    loadedColorsFromTxt = true;
                }
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

        params.centralMass = 1.0;      // 1 масса Солнца
        params.baseMass    = 3.0e-6;   // Масса Земли в солнечных
        params.minRadius   = 5.0f;     // AU
        params.maxRadius   = 40.0f;    // AU
        params.spread      = 2.0f;     // AU

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
            // В питоне ComboBox index начинается с 0, а у нас в консоли с 1. Переводим
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
    
    BodySystem bodySystem(objs); // создание сохрянем информацию о системе
    bodySystem.transPointToSystem(objs); // переходим в систему объектов

    const double initialEnergy = physics::calculateTotalEnergy(objs);
    double energyCache = initialEnergy;

    Camera cam;
    SimState state; 

    Control control(renderer.getWindow(), objs, cam, state);
    control.attach();
    int counter = 0;
    double lastTime = glfwGetTime();
    double accumulator = 0.0;

    bool isRealTime;
    if (!is_cli) {
        std::cout << "Режим реального времени? [1 - Да / 0 Нет]: " << std::flush;
        std::cin >> isRealTime;
    } else {
        isRealTime = (getCmdOption(argc, argv, "--realtime", "1") == "1");
    }

    int stepCounter = 0;
    int frameCounter = 0;

    if (isRealTime) {
        while (!glfwWindowShouldClose(renderer.getWindow()) && state.running && !g_Interrupt) {
            double now = glfwGetTime();
            double frameRealDt = now - lastTime;
            
            state.deltaTime = frameRealDt;
            lastTime = now;
            frameRealDt *= state.timeScale;
            if (frameRealDt > 0.1) frameRealDt = 0.1;
            if (frameRealDt < 0.0) frameRealDt = 0.0;
            // Зажатый Shift: время симуляции не идёт (удобно летать камерой / менять Shift+− скорость камеры).
            GLFWwindow* win = renderer.getWindow();
            const bool shiftFreeze = (glfwGetKey(win, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
                || (glfwGetKey(win, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS);
            if (!shiftFreeze) {
                accumulator += frameRealDt;
            }

            int substeps = 0;

            constexpr int MAX_SUBSTEPS = 100000;
            while (accumulator >= fixedDt && substeps < MAX_SUBSTEPS && !g_Interrupt) {
                simulationStep(objs, fixedDt, state.pause, false);
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

            control.updateCameraFromKeys();

            renderer.renderFrame(objs, graphics, cam);
            ++frameCounter;

            double dtForFps = frameRealDt / std::max(1.0f, state.timeScale);
            double fps = (dtForFps > 0.0) ? 1.0 / dtForFps : 0.0;
            
            std::snprintf(title, sizeof(title),
                  "REAL-TIME | Sim: %.1fx | Cam: %.2fx | FPS: %.0f | Obj: %zu | Time: %.2f%s",
                  state.timeScale, control.getCameraMoveScale(), fps, objs.size(), gSimTime,
                  shiftFreeze ? " | [SHIFT freeze]" : "");
            glfwSetWindowTitle(renderer.getWindow(), title);

            glfwPollEvents();
        }

    } else {
        int startTime = static_cast<int>(time(NULL));
        double targetTime = 0.0;
        
        {
            const double yearsBetweenSaves = fixedDt * static_cast<double>(saveIntervalSteps);
            const double mbPerSimYear =
                (yearsBetweenSaves > 0.0 && !objs.empty())
                    ? (static_cast<double>(objs.size()) * 3.0 * static_cast<double>(sizeof(float))
                       / (yearsBetweenSaves * 1024.0 * 1024.0))
                    : 0.0;
                    
            if (!is_cli) {
                std::cout << std::fixed << std::setprecision(6);
                std::cout << "Просчёт (не в реальном времени). Один шаг dt = " << fixedDt << " года.\n";
                std::cout << "Между сохранёнными кадрами: " << yearsBetweenSaves << " года ("
                          << saveIntervalSteps << " шаг(ов)). Ориентир размера треков: ~" << std::setprecision(3)
                          << mbPerSimYear << " МБ на 1 год симуляции, ~" << std::setprecision(6)
                          << " (оценка по float-позициям; deflate в HDF5 обычно меньше).\n";
                std::cout << std::setprecision(15);
                std::cout << "На какое симуляционное время просчитываем? (годы): " << std::flush;
                std::cin >> targetTime;
            } else {
                targetTime = std::stod(getCmdOption(argc, argv, "--target", "10.0"));
            }
            
            if (targetTime > 0.0 && mbPerSimYear > 0.0 && !is_cli) {
                std::cout << std::fixed << std::setprecision(2)
                          << "Ориентир объёма треков на весь прогон: ~" << (mbPerSimYear * targetTime)
                          << " МБ (грубо, без учёта initial_data и сжатия).\n";
            }
        }
        
        double bakeStartRealTime = glfwGetTime();

        glfwSwapInterval(0);
        std::cout << "Начинаем расчет..." << std::endl;

        while (!glfwWindowShouldClose(renderer.getWindow()) && state.running && gSimTime < targetTime && !g_Interrupt) {
            simulationStep(objs, fixedDt, false, false);
            gSimTime += fixedDt;
            stepCounter += 1;
            
            if (stepCounter % (1 * 10) == 0) {
                glfwPollEvents();

                double progress = (gSimTime / targetTime) * 100.0;
                double elapsed = glfwGetTime() - bakeStartRealTime;
                double progress_dec = gSimTime / targetTime;
                double etaSeconds = 0.0;
                if (progress_dec > 0.001) {
                    etaSeconds = (elapsed / progress_dec) - elapsed;
                    if (etaSeconds < 0.0) etaSeconds = 0.0;
                }
                int etaTotalSec = static_cast<int>(etaSeconds);
                int etaMin = etaTotalSec / 60;
                int etaSec = etaTotalSec % 60;
                
                std::snprintf(title, sizeof(title),
                      "BAKING... %.1f%% | ETA: %02d:%02d | Time: %.2f / %.2f | Saved: %zu",
                      progress, etaMin, etaSec, gSimTime, targetTime, frameIndex);
                glfwSetWindowTitle(renderer.getWindow(), title);
                std::cout << "PROGRESS:" << progress << std::endl;
            }
            if (stepCounter % saveIntervalSteps == 0) {
                WriteSimulationFrame(framesFile, objs, graphics, frameIndex);
                ++frameIndex;
            }
        }
        std::cout << "Расчет завершен за " << (static_cast<int>(time(NULL)) - startTime) << " сек.\n";
    }

    if (g_Interrupt) {
        std::cout << "Получен SIGINT (Ctrl+C). Завершаем симуляцию и сохраняем текущее состояние...\n";
    }

    simulationStep(objs, 0.0, true, true);
    const std::string finalTextPath = "data/" + filename + "_final.txt";
    if (SaveSystemToTextFile(finalTextPath, objs)) {
        std::cout << "Финальное состояние сохранено в " << finalTextPath << "\n";
    } else {
        std::cout << "Не удалось сохранить финальное состояние в " << finalTextPath << "\n";
    }

    CloseSimulationFile(framesFile, frameIndex);
    return 0;
}