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
#include <random>
#include <filesystem>
#include "constants.hpp"
#include "object.hpp"
#include "renderer.hpp"
#include "control.hpp"
#include "bodysystem.hpp"
#include "../data/data.hpp"
#include "barnes_hut.hpp"

#ifndef VEL_SCALE
#define VEL_SCALE 1.0f
#endif

#ifdef _WIN32
#include <windows.h>
#endif

struct AppState {
    bool running   = true;
    bool pause     = false;

    glm::vec3 cameraPos   = glm::vec3(0.0f, 50.0f, 250.0f);
    glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
    glm::vec3 cameraUp    = glm::vec3(0.0f, 1.0f, 0.0f);

    float lastX = WINDOW_W / 2.0f; 
    float lastY = WINDOW_H / 2.0f;
    float yaw   = -90.0f;
    float pitch = 0.0f;

    float dt = 0.0f;
    float timeScale = 1.0f; 
    
    const float fixedDt = FIXED_DT; 
    float initMass = 5.0f * std::pow(10.0f, 20.0f); 
    
    std::vector<Object> objs = {};
    char title[128];
};

static inline float calculateRadius(double m) {
    const double r_m = std::cbrt((3.0 * m) / (4.0 * glm::pi<double>() * RHO_DEFAULT));
    return static_cast<float>(r_m / SCALE_FACTOR); 
}

static inline float RadiusKm(double m, double rho) {
    const double r_m = cbrt((3.0 * m) / (4.0 * 3.14159265358979323846 * rho));
    return static_cast<float>(r_m / 50000.0); 
}

void colorFromMass(std::vector<Object>& objs) {
    double minMass = FLT_MAX;
    double maxMass = -FLT_MAX;

    for (const auto& obj : objs) {
        minMass = std::min(minMass, obj.getMass());
        maxMass = std::max(maxMass, obj.getMass());
    }

    const glm::vec3 burgundy = glm::vec3(0.50f, 0.00f, 0.125f);
    const glm::vec3 white    = glm::vec3(1.00f, 1.00f, 1.00f);

    double denom = std::log10(std::max(1e-30, maxMass)) - std::log10(std::max(1e-30, minMass));
    if (!std::isfinite(denom) || std::abs(denom) < 1e-12) {
        for (auto& o : objs) o.setColor(glm::vec4(glm::mix(burgundy, white, 0.5f), 1.0f)); 
        return;
    }

    for (auto& obj : objs) {
        float m = glm::max(obj.getMass(), 1e-30); 
        float t = (std::log10(m) - std::log10(minMass)) / denom;
        t = glm::clamp(t, 0.0f, 1.0f);
        glm::vec3 rgb = glm::mix(burgundy, white, t);
        obj.setColor(glm::vec4(rgb, 1.0f)); 
    }
}

void spawnSystem(std::vector<Object>& out, int N, double centralMass, double satMassBase,
                 float rMin_km, float rMax_km, unsigned seed = 42)
{
    out.clear();
    out.reserve(static_cast<size_t>(N) + 1);
    out.emplace_back(glm::dvec3(0.0), glm::dvec3(0.0), centralMass, glm::dvec3(0.0));
    out.back().setColor(glm::vec4(1.0f, 1.0f, 1.0f, 1.0f)); 
    if (N <= 0) return;
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> u01(0.0f, 1.0f);
    std::uniform_real_distribution<float> uAngle(0.0f, glm::pi<float>() * 2.0f);
    const float ringWidth = std::max(0.0f, rMax_km - rMin_km);
    const float rJitterMax = 0.01f * ringWidth;
    for (int i = 0; i < N; ++i) {
        float t = (static_cast<float>(i) + 0.5f) / static_cast<float>(N);
        float r2_min = rMin_km * rMin_km;
        float r2_max = rMax_km * rMax_km;
        float r = std::sqrt(r2_min + (r2_max-r2_min)*t);
        r += (u01(rng) * 2.0f - 1.0f) * rJitterMax;
        r = std::clamp(r, rMin_km, rMax_km);
        float a = uAngle(rng);
        glm::dvec3 pos(r * std::cos(a), 0.0, r * std::sin(a)); 
        glm::dvec3 tdir(-std::sin(a), 0.0, std::cos(a)); 
        double v_circ_mps = std::sqrt((G * centralMass) / (static_cast<double>(r) * 1000.0));
        double v_kmps     = v_circ_mps / 1000.0 * VEL_SCALE;
        out.emplace_back(pos, tdir * v_kmps, satMassBase, glm::dvec3(0.0));
    }
    colorFromMass(out);
}

void simulationStepBrutForceCPU(std::vector<Object>& objs, float dt, bool pause){
    
    for (auto& obj : objs) {
        obj.setAcc(glm::dvec3(0.0));
    }
    
    // 2. Расчет сил
    for (size_t i = 0; i < objs.size(); ++i) {
        Object& obj = objs[i];
        for (size_t j = i + 1; j < objs.size(); ++j) {
            Object& obj2 = objs[j];

            glm::dvec3 delta = obj2.getPos() - obj.getPos(); 
            double distance = glm::length(delta); 

            if (distance <= 0.0) continue;

            glm::dvec3 dir = delta / distance;
            
            double radius1 = calculateRadius(obj.getMass());
            double radius2 = calculateRadius(obj2.getMass());
            double combinedRadius = radius1 + radius2;

            double effectiveDistance = std::max(distance, combinedRadius * 0.1); 
            double dist_m = effectiveDistance * 1000.0;        
            
            double F = (G * obj.getMass() * obj2.getMass()) / (dist_m * dist_m); 
            
            double acc1_kmps2 = (F / obj.getMass()) / 1000.0; 
            double acc2_kmps2 = (F / obj2.getMass()) / 1000.0;
            
            obj.setAcc(obj.getAcc() + dir * acc1_kmps2);
            obj2.setAcc(obj2.getAcc() - dir * acc2_kmps2);
        }
    }

    for (size_t i = 0; i < objs.size(); ++i) {
        Object& obj = objs[i];
        
        if (!pause) {
            glm::dvec3 newVel = obj.getVel() + obj.getAcc() * static_cast<double>(dt);
            obj.setVel(newVel);
            
            glm::dvec3 newPos = obj.getPos() + newVel * static_cast<double>(dt);
            obj.setPos(newPos);
        }
        
        for (size_t j = i + 1; j < objs.size(); ++j) {
            Object& obj2 = objs[j];
            
            glm::dvec3 delta = obj2.getPos() - obj.getPos();
            double distance = glm::length(delta);
            
            double radius1 = calculateRadius(obj.getMass());
            double radius2 = calculateRadius(obj2.getMass());
            double combinedRadius = radius1 + radius2;

            if (distance < combinedRadius) {
                glm::dvec3 normal = delta / distance;
                glm::dvec3 relativeVelocity = obj.getVel() - obj2.getVel(); 
                double relVelAlongNormal = glm::dot(relativeVelocity, normal);

                if (relVelAlongNormal < 0.0) {
                    double restitution = 0.8;
                    double invMass1 = 1.0 / obj.getMass(); 
                    double invMass2 = 1.0 / obj2.getMass();
                    double impulseScalar = -(1.0 + restitution) * relVelAlongNormal / (invMass1 + invMass2);
                    glm::dvec3 impulse = normal * impulseScalar;
                    
                    obj.setVel(obj.getVel() + impulse * invMass1);
                    obj2.setVel(obj2.getVel() - impulse * invMass2);
                }

                double penetration = combinedRadius - distance;
                if (penetration > 0.0) {
                    double invMass1 = 1.0 / obj.getMass();
                    double invMass2 = 1.0 / obj2.getMass();
                    double invMassSum = invMass1 + invMass2;
                    if (invMassSum > 0.0) {
                        double correctionScale = penetration / invMassSum;
                        glm::dvec3 correction = normal * correctionScale;
                        
                        obj.setPos(obj.getPos() - correction * invMass1);
                        obj2.setPos(obj2.getPos() + correction * invMass2);
                    }
                }
            }
        }
    }
}


GLFWwindow* StartGLU();


int main() {

#ifdef _WIN32
    SetConsoleOutputCP(65001);
    SetConsoleCP(65001);
#endif
    GLFWwindow* window = StartGLU();
    if (!window) {
        std::cerr << "Window or OpenGL context creation failed.\n";
        return -1;
    }
    Renderer renderer(WINDOW_W, WINDOW_H); 
    renderer.setProjection(FOV, (float)WINDOW_W/(float)WINDOW_H, Z_NEAR, Z_FAR);
    using Handler = void(*)(std::vector<Object>& objs, float dt, bool pause);
    Handler simulationStep = nullptr;
    cameraPos = glm::vec3(0.0f, 50.0f, 250.0f);
    bool loaded = false;
    char mode;
    std::cout << "Выберите алгоритм: брутфорс или Барнс-Хат? [0/1]: " << std::flush;
    std::cin >> mode;
    if (mode == '0') {
        simulationStep = &simulationStepBrutForceCPU;
    }
    else{
        simulationStep = &simulationStepBarnesHutCPU;
    }
    std::cout << "Загружаем сценарий из HDF5 или генерируем систему рандомно? [0/1]: " << std::flush;
    std::cin >> mode;
    if (mode == '0'){
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
                colorFromMass(objs);
                loaded = true;
                std::cout << "Loaded " << objs.size() << " objects\n";
            }

        }
    }
    if (!loaded){
        double M_central  = static_cast<double>(initMass) * 1000;
        double M_sat_base = static_cast<double>(initMass);
        int numObjs;
        std::cout << "Сколько тел загружаем?: " << std::flush;
        std::cin >> numObjs;
        spawnSystem(objs, numObjs, M_central, M_sat_base, /*rMin*/300.0f, /*rMax*/70000.0f, /*seed*/42);
    }
    // Чтение с HDF5
    BodySystem bodySystem(objs); // создание сохрянем информацию о системе
    bodySystem.transPointToSystem(objs); // переходим в систему объектов
    // Управление
    Control control(window, state.objs,
                    state.cameraPos, state.cameraFront, state.cameraUp,
                    state.dt, state.timeScale, state.pause, state.running,
                    state.yaw, state.pitch, state.lastX, state.lastY,
                    state.initMass);
    control.attach();

    double lastTime = glfwGetTime();
    double accumulator = 0.0;
    while (!glfwWindowShouldClose(window) && state.running) {
        double now = glfwGetTime();
        double frameRealDt = now - lastTime;
        state.dt = (float)frameRealDt; 
        lastTime = now;
        frameRealDt *= state.timeScale;
        accumulator += frameRealDt;
        int substeps = 0;
        const int MAX_SUBSTEPS = 8;
        while (accumulator >= state.fixedDt && substeps < MAX_SUBSTEPS) {
            simulationStep(state.objs, state.fixedDt, state.pause);
            accumulator -= state.fixedDt;
            ++substeps;
        }
        while (accumulator >= state.fixedDt) {
            simulationStep(state.objs, state.fixedDt, state.pause);
            accumulator -= state.fixedDt;
        }
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        renderer.updateView(state.cameraPos, state.cameraFront, state.cameraUp);
        renderer.drawObjects(state.objs);
        double dtForFps = frameRealDt / std::max(1.0f, state.timeScale);
        double fps = (dtForFps > 0.0) ? 1.0 / dtForFps : 0.0;
        std::snprintf(state.title, sizeof(state.title),
              "3D_TEST | timeScale: %.2fx | FPS: %.0f | Objects =: %zu | Paused: %s",
              state.timeScale, fps, state.objs.size(), state.pause ? "YES" : "NO");
        glfwSetWindowTitle(window, state.title);
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    glfwTerminate();
    return 0;
}

GLFWwindow* StartGLU() {
    if (!glfwInit()) {
        std::cout << "Failed to initialize GLFW\n";
        return nullptr;
    }
    glfwInitHint(GLFW_PLATFORM, GLFW_PLATFORM_X11);
    GLFWwindow* window = glfwCreateWindow(WINDOW_W, WINDOW_H, "3D_TEST", nullptr, nullptr);
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
    glViewport(0, 0, WINDOW_W, WINDOW_H);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    return window;
}