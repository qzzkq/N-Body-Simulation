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
#include "object.hpp"
#include "renderer.hpp"
#include "control.hpp"
#include "bodysystem.hpp"
#include "data.hpp"
#include "barnes_hut.hpp"

#ifndef VEL_SCALE
#define VEL_SCALE 1.0f
#endif

#ifdef _WIN32
#include <windows.h>
#endif


bool running = true;
bool pause   = false;

glm::vec3 cameraPos   = glm::vec3(0.0f, 0.0f, 1.0f);
glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
glm::vec3 cameraUp    = glm::vec3(0.0f, 1.0f, 0.0f);

float lastX = 400.0f, lastY = 300.0f;
float yaw   = -90.0f;
float pitch = 0.0f;

float dt = 0.0f;
float lastFrame = 0.0f;

float timeScale = 1.0f; // переменная для ускорения/замедления времени
float fixedDt = 1.0f / 240.0f; // шаг времени;
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

static inline float RadiusKm(double m, double rho) {
    const double r_m = cbrt((3.0 * m) / (4.0 * 3.14159265358979323846 * rho));
    return static_cast<float>(r_m / 50000.0); // ← было /100000.0
}

void colorFromMass(std::vector<Object>& objs) {
    double minMass = FLT_MAX;
    double maxMass = -1;

    for (size_t i = 0; i < objs.size() ; i++){
        Object& obj = objs[i];
        minMass = (minMass > obj.mass) ? obj.mass : minMass;
        maxMass = (maxMass < obj.mass) ? obj.mass : maxMass;
    }

    const glm::vec3 burgundy = glm::vec3(0.50f, 0.00f, 0.125f);
    const glm::vec3 white    = glm::vec3(1.00f, 1.00f, 1.00f);

    double denom = std::log10(std::max(1e-30, minMass)) - std::log10(std::max(1e-30, maxMass));
    if (!std::isfinite(denom) || std::abs(denom) < 1e-12) {
        for (auto& o : objs) o.color = glm::vec4(glm::mix(burgundy, white, 0.5f), 1.0f);
        return;
    }

    for (size_t i = 0; i < objs.size(); i++){
        Object& obj = objs[i];
        float m = glm::max(obj.mass, 1e-30);
        float t = (std::log10(m) - std::log10(minMass)) /
                (std::log10(maxMass) - std::log10(minMass));
        t = glm::clamp(t, 0.0f, 1.0f);
        glm::vec3 rgb = glm::mix(burgundy, white, t);
        obj.color = glm::vec4(rgb, 1.0f);
    }
}

void spawnSystem(std::vector<Object>& out, int N, double centralMass, double satMassBase,
                 float rMin_km, float rMax_km, unsigned seed /*= 42*/)
{
    out.clear();
    out.reserve(static_cast<size_t>(N) + 1);

    Object center(glm::vec3(0), glm::vec3(0), centralMass, 141000.0f, std::nullopt);

    center.Initalizing = false;
    center.radius = RadiusKm(center.mass, center.density);
    center.UpdateVertices();
    out.push_back(center);

    if (N <= 0) return;

    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> u01(0.0f, 1.0f);
    std::uniform_real_distribution<float> uAngle(0.0f, 6.28318530718f);

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

        glm::vec3 pos(r * std::cos(a), 0.0f, r * std::sin(a));
        glm::vec3 tdir(-std::sin(a), 0.0f, std::cos(a));

        double v_circ_mps = std::sqrt((G * centralMass) / (static_cast<double>(r) * 1000.0));
        float  v_kmps     = static_cast<float>(v_circ_mps / 1000.0f) * VEL_SCALE;

        Object o(pos, tdir * v_kmps, satMassBase, 1410.0f, std::nullopt);
        o.Initalizing = false;
        o.radius = RadiusKm(o.mass, o.density);
        o.UpdateVertices();

        out.push_back(std::move(o));
    }

    colorFromMass(out);

}


void simulationStepBrutForceCPU(std::vector<Object>& objs, float dt, bool pause){
    for (size_t i = 0; i < objs.size(); ++i) {
        Object& obj = objs[i];
        if (obj.Initalizing) {
            obj.radius = std::pow((3.0 * obj.mass / obj.density) / (4.0 * 3.14159265359), 1.0/3.0) / 100000.0;
            obj.UpdateVertices();
        }

        for (size_t j = i + 1; j < objs.size(); ++j) {
            Object& obj2 = objs[j];
            if (obj.Initalizing || obj2.Initalizing) {
                continue;
            }

            glm::dvec3 delta = obj2.GetPos() - obj.GetPos();
            double distance = std::sqrt(delta.x * delta.x + delta.y * delta.y + delta.z * delta.z);
            if (distance <= 0.0) {
                continue;
            }

            glm::dvec3 dir = delta / distance;
            double combinedRadius = obj.radius + obj2.radius;

            double effectiveDistance = std::max(distance, combinedRadius);
            double dist_m = effectiveDistance * 1000.0;        // м
            double F = (G * obj.mass * obj2.mass) / (dist_m * dist_m);              // Н
            float acc1_kmps2 = static_cast<float>((F / obj.mass)  / 1000.0);        // км/с²
            float acc2_kmps2 = static_cast<float>((F / obj2.mass) / 1000.0);        // км/с²
            glm::vec3 accObj  = glm::vec3(dir * static_cast<double>(acc1_kmps2));
            glm::vec3 accObj2 = glm::vec3(-dir * static_cast<double>(acc2_kmps2));

            if (!pause) {
                obj.accelerate(accObj.x, accObj.y, accObj.z, dt);
                obj2.accelerate(accObj2.x, accObj2.y, accObj2.z, dt);
            }

            if (distance < combinedRadius) {
                glm::vec3 normal = glm::vec3(dir);
                glm::vec3 relativeVelocity = glm::vec3(obj.velocity - obj2.velocity);
                float relVelAlongNormal = glm::dot(relativeVelocity, normal);

                if (relVelAlongNormal < 0.0f) {
                    double restitution = 0.8;
                    double invMass1 = 1.0 / obj.mass;
                    double invMass2 = 1.0 / obj2.mass;
                    double impulseScalar = -(1.0 + restitution) * static_cast<double>(relVelAlongNormal) / (invMass1 + invMass2);
                    glm::dvec3 impulse = glm::dvec3(normal) * impulseScalar;
                    obj.velocity += impulse * invMass1;
                    obj2.velocity -= impulse * invMass2;
                }

                double penetration = combinedRadius - distance;
                if (penetration > 0.0) {
                    double invMass1 = 1.0 / obj.mass;
                    double invMass2 = 1.0 / obj2.mass;
                    double invMassSum = invMass1 + invMass2;
                    if (invMassSum > 0.0) {
                        double correctionScale = penetration / invMassSum;
                        glm::dvec3 correction = glm::dvec3(normal) * correctionScale;
                        obj.position -= correction * invMass1;
                        obj2.position += correction * invMass2;
                    }
                }
            }
        }

        if (!pause) {
            obj.UpdatePos(dt);
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

    Renderer renderer(800, 600, vertexShaderSource, fragmentShaderSource);
    renderer.setProjection(65.0f, 800.0f/600.0f, 8.3f, 100000.0f);
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
    Control control(window, objs,
                    cameraPos, cameraFront, cameraUp,
                    dt, timeScale, pause, running,
                    yaw, pitch, lastX, lastY,
                    initMass);
    control.attach();

    double lastTime = glfwGetTime();
    double accumulator = 0.0;
    while (!glfwWindowShouldClose(window) && running) {
        double now = glfwGetTime();
        double frameRealDt = now - lastTime;
        dt = frameRealDt;
        lastTime = now;
        frameRealDt *= timeScale;
        accumulator += frameRealDt;

        int substeps = 0;
        const int MAX_SUBSTEPS = 8;
        while (accumulator >= fixedDt && substeps < MAX_SUBSTEPS) {
            simulationStep(objs, fixedDt, pause);
            accumulator -= fixedDt;
            ++substeps;
        }
        while (accumulator >= fixedDt) {
            simulationStep(objs, fixedDt, pause);
            accumulator -= fixedDt;
        }

        // Отрисовка
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        renderer.updateView(cameraPos, cameraFront, cameraUp);
        renderer.drawObjects(objs);

        double dtForFps = frameRealDt / std::max(1.0f, timeScale);
        double fps = (dtForFps > 0.0) ? 1.0 / dtForFps : 0.0;
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
    if (!glfwInit()) {
        std::cout << "Failed to initialize GLFW\n";
        return nullptr;
    }
    glfwInitHint(GLFW_PLATFORM, GLFW_PLATFORM_X11);
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