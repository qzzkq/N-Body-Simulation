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
#include "physics.hpp"
using namespace H5;

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

struct InitRecord {
    glm::dvec3 position;
    double     mass;
    double     radius;
};

static CompType GetInitType() {
    hsize_t v3dims[1] = {3};
    ArrayType vec3Type(PredType::NATIVE_DOUBLE, 1, v3dims);
    CompType t(sizeof(InitRecord));
    t.insertMember("position", HOFFSET(InitRecord, position), vec3Type);
    t.insertMember("mass",     HOFFSET(InitRecord, mass),     PredType::NATIVE_DOUBLE);
    t.insertMember("radius",   HOFFSET(InitRecord, radius),   PredType::NATIVE_DOUBLE);
    return t;
}

static ArrayType GetPosType() {
    hsize_t v3dims[1] = {3};
    return ArrayType(PredType::NATIVE_DOUBLE, 1, v3dims);
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

    H5File file;
    try {
        file.openFile(fileName, H5F_ACC_RDONLY);
    } catch(...) {
        std::cerr << "Failed to open file.\n";
        return 1;
    }

    hsize_t numBodiesH5 = 0;
    hsize_t numFramesH5 = 0;
    double fileDt = 1.0/10;

    try {
        file.openAttribute("num_bodies").read(PredType::NATIVE_HSIZE, &numBodiesH5);
        file.openAttribute("num_frames").read(PredType::NATIVE_HSIZE, &numFramesH5);
        file.openAttribute("dt").read(PredType::NATIVE_DOUBLE, &fileDt);
    } catch (...) {
        std::cerr << "Error File header missing.\n";
        return 1;
    }

    size_t numBodies = (size_t)numBodiesH5;
    size_t numFrames = (size_t)numFramesH5;

    std::cout << "Opened: " << fileName << "\n"
              << "  Bodies: " << numBodies << "\n"
              << "  Frames: " << numFrames << "\n"
              << "  DT: " << fileDt << " s\n";

    if (numBodies == 0 || numFrames == 0) return 1;

    objs.clear();
    objs.reserve(numBodies);
    try {
        DataSet initSet = file.openDataSet("initial_data");
        std::vector<InitRecord> initBuf(numBodies);
        initSet.read(initBuf.data(), GetInitType());

        for (const auto& rec : initBuf) {
            Object o(rec.position, glm::dvec3(0), rec.mass, 1410.0f, std::nullopt);
            o.Initalizing = false;
            o.radius = (float)rec.radius;
            objs.push_back(o);
        }
    } catch(...) {
        std::cerr << "Error reading /initial_data.\n";
        return 1;
    }
    physics::colorFromMass(objs);

    DataSet tracksSet;
    try {
        tracksSet = file.openDataSet("tracks");
    } catch(...) {
        std::cerr << "Error opening /tracks.\n";
        return 1;
    }

    GLFWwindow* window = StartGLU();
    if (!window) return 1;

    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    Renderer renderer(width, height, vertexShaderSource, fragmentShaderSource);
    renderer.setProjection(65.0f, (float)width/height, 0.1f, 100000.0f);

    cameraPos   = glm::vec3(0.0f, 50.0f, 250.0f);
    cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
    cameraUp    = glm::vec3(0.0f, 1.0f,  0.0f);

    renderer.setRenderMode(RenderMode::Sphere);
    
    float initMass = 1.0f;

    Control control(window, objs, cameraPos, cameraFront, cameraUp, dt, timeScale, pause, running, yaw, pitch, lastX, lastY, initMass);
    control.attach();

    double playbackTime = 0.0;
    double maxTime = (numFrames - 1) * fileDt;

    double lastRealTime = glfwGetTime();

    size_t currentFrameIdx = 0;
    size_t lastReadFrameIdx = 99999999999;
    size_t pageNumber = 0;
    
    std::vector<glm::dvec3> posBuffer(numBodies);

    while (!glfwWindowShouldClose(window) && running) {
        double now = glfwGetTime();
        double realDt = now - lastRealTime;
        lastRealTime = now;
        dt = (float)realDt;

        if (!pause) {
            playbackTime += realDt * timeScale;
        }

        if (playbackTime > maxTime) playbackTime = 0.0;
        if (playbackTime < 0.0) playbackTime = maxTime;

        currentFrameIdx = (size_t)(playbackTime / fileDt);
        if (currentFrameIdx >= numFrames) currentFrameIdx = numFrames - 1;

        if (currentFrameIdx != lastReadFrameIdx) {
            hsize_t offset[2] = { static_cast<hsize_t>(currentFrameIdx), 0 };
            hsize_t count[2]  = { 1, static_cast<hsize_t>(numBodies) };
            
            DataSpace fileSpace = tracksSet.getSpace();
            fileSpace.selectHyperslab(H5S_SELECT_SET, count, offset);

            hsize_t memDims[2] = { 1, static_cast<hsize_t>(numBodies) };
            DataSpace memSpace(2, memDims);

            tracksSet.read(posBuffer.data(), GetPosType(), memSpace, fileSpace);

            for (size_t i = 0; i < numBodies; ++i) {
                objs[i].position = posBuffer[i];
            }

            lastReadFrameIdx = currentFrameIdx;
        }

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        renderer.updateView(cameraPos, cameraFront, cameraUp);
        renderer.drawObjects(objs);

        std::snprintf(title, sizeof(title),
                      "Replay | Time: %.2f / %.2f | Frame: %zu | Speed: x%.2f",
                      playbackTime, maxTime, currentFrameIdx, timeScale);
        glfwSetWindowTitle(window, title);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}
