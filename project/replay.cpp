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
#include "state.hpp"
#include "camera.hpp"
using namespace H5;

char title[256];

std::vector<Object> objs;

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


const int INTERPOLATION_STEPS = 10;
int main(int argc, char** argv) {

    Camera cam; 
    SimState state; 

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
    double fixedDt = 1.0f/60;

    try {
        file.openAttribute("num_bodies").read(PredType::NATIVE_HSIZE, &numBodiesH5);
        file.openAttribute("num_frames").read(PredType::NATIVE_HSIZE, &numFramesH5);
        file.openAttribute("dt").read(PredType::NATIVE_DOUBLE, &fixedDt);
    } catch (...) {
        std::cerr << "Error File header missing.\n";
        return 1;
    }

    size_t numBodies = (size_t)numBodiesH5;
    size_t numFrames = (size_t)numFramesH5;

    std::cout << "Opened: " << fileName << "\n"
              << "  Bodies: " << numBodies << "\n"
              << "  Frames: " << numFrames << "\n"
              << "  DT: " << fixedDt << " s\n";

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

    bool fullscreen = false;
    bool maximized = true;
    Renderer renderer;
    if (!renderer.init(1280, 720, "N-Body simulation", fullscreen, maximized)) {
        return 1;
    }

    renderer.setProjection(65.0f, 0.1f, 100000.0f);

    cam.pos = glm::vec3(0.0f, 50.0f, 250.0f);
    cam.front = glm::vec3(0.0f, 0.0f, -1.0f);
    cam.up    = glm::vec3(0.0f, 1.0f,  0.0f);

    renderer.setRenderMode(RenderMode::Sphere);

    Control control(renderer.getWindow(), objs, cam, state);
    control.attach();

    double playbackTime = 0.0;
    double maxTime = (numFrames - 1) * fixedDt; 

    double lastRealTime = glfwGetTime();

    size_t frameA = 0;
    size_t lastReadFrameIdx = 99999999999;
    size_t pageNumber = 0;
    
    std::vector<glm::dvec3> posBufferA(numBodies);
    std::vector<glm::dvec3> posBufferB(numBodies);
    std::vector<glm::dvec3> targetPosition(numBodies);
    double ratio;
    while (!glfwWindowShouldClose(renderer.getWindow()) && state.running) {
        double now = glfwGetTime();
        double realDt = now - lastRealTime;
        lastRealTime = now;
        state.deltaTime = (float)realDt;

        if (!state.pause) {
            playbackTime += realDt * state.timeScale;
        }

        if (playbackTime > maxTime) playbackTime = 0.0;
        if (playbackTime < 0.0) playbackTime = maxTime;
        ratio = playbackTime / fixedDt;
        frameA = (size_t)ratio;
        size_t frameB = frameA + 1;
        if (frameA >= numFrames) frameA = numFrames - 1;
        if (frameB >= numFrames) frameB = numFrames - 1;

        if (frameA < lastReadFrameIdx && lastReadFrameIdx != 99999999999) {
            for (auto& obj : objs) {
                obj.trail.clear(); 
            }
        }

        if (lastReadFrameIdx != frameA) {
            hsize_t count[2] = { 1, static_cast<hsize_t>(numBodies) };
            hsize_t memDims[2] = { 1, static_cast<hsize_t>(numBodies) };
            DataSpace memSpace(2, memDims);

            hsize_t offsetA[2] = { static_cast<hsize_t>(frameA), 0 };
            DataSpace fileSpaceA = tracksSet.getSpace();
            fileSpaceA.selectHyperslab(H5S_SELECT_SET, count, offsetA);
            tracksSet.read(posBufferA.data(), GetPosType(), memSpace, fileSpaceA);

            hsize_t offsetB[2] = { static_cast<hsize_t>(frameB), 0 };
            DataSpace fileSpaceB = tracksSet.getSpace();
            fileSpaceB.selectHyperslab(H5S_SELECT_SET, count, offsetB);
            tracksSet.read(posBufferB.data(), GetPosType(), memSpace, fileSpaceB);
        }

        //std::vector<glm::dvec3> targetPosition = glm::mix(posBufferA, posBufferB, (float)(ratio - frameA) / (float)(frameB - frameA));

        for (size_t i = 0; i < numBodies; ++i) {
            targetPosition[i] = glm::mix(posBufferA[i], posBufferB[i], (float)(ratio - frameA));
        }

        for (size_t i = 0; i < numBodies; ++i) {
            objs[i].position = targetPosition[i];
        }

        // for(auto& obj : objs) obj.updateTrail(); 

        lastReadFrameIdx = frameA;

    renderer.renderFrame(objs, cam); 

        std::snprintf(title, sizeof(title),
                      "Replay | Time: %.2f / %.2f | Frame: %zu | Speed: x%.2f",
                      playbackTime, maxTime, frameA, state.timeScale);
        glfwSetWindowTitle(renderer.getWindow(), title);
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
    
}
