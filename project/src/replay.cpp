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
#include <utility>
#include <array>

#ifdef _WIN32
#include <windows.h>
#endif

#include "object.hpp"
#include "graphic_state.hpp"
#include "renderer.hpp"
#include "control.hpp"
#include "physics.hpp"
#include "state.hpp"
#include "camera.hpp"
using namespace H5;

char title[256];

std::vector<Object> objs;
std::vector<GraphicState> graphics;

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

struct PosRecord {
    glm::dvec3 position;
};

struct PosRecordFile {
    glm::vec3 position;
};

static CompType GetTrackType() {
    hsize_t v3dims[1] = {3};
    ArrayType vec3Type(PredType::NATIVE_FLOAT, 1, v3dims);
    CompType t(sizeof(PosRecordFile));
    t.insertMember("position", HOFFSET(PosRecordFile, position), vec3Type);
    return t;
}

static CompType GetTrackMemType() {
    hsize_t v3dims[1] = {3};
    ArrayType vec3Type(PredType::NATIVE_DOUBLE, 1, v3dims);
    CompType t(sizeof(PosRecord));
    t.insertMember("position", HOFFSET(PosRecord, position), vec3Type);
    return t;
}

static glm::dvec3 catmullRom(const glm::dvec3& p0,
                             const glm::dvec3& p1,
                             const glm::dvec3& p2,
                             const glm::dvec3& p3,
                             double t) {
    const double t2 = t * t;
    const double t3 = t2 * t;
    return 0.5 * ((2.0 * p1) +
                  (-p0 + p2) * t +
                  (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t2 +
                  (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t3);
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
    graphics.resize(objs.size());
    bool loadedReplayColors = false;
    try {
        DataSet colorSet = file.openDataSet("body_colors");
        DataSpace csp = colorSet.getSpace();
        int crank = csp.getSimpleExtentNdims();
        if (crank == 2) {
            hsize_t cdims[2] = {0, 0};
            csp.getSimpleExtentDims(cdims);
            if (cdims[0] == static_cast<hsize_t>(numBodies) && cdims[1] == 4) {
                std::vector<float> cbuf(numBodies * 4);
                colorSet.read(cbuf.data(), PredType::NATIVE_FLOAT);
                for (size_t i = 0; i < numBodies; ++i) {
                    graphics[i].color = glm::vec4(
                        cbuf[i * 4 + 0], cbuf[i * 4 + 1], cbuf[i * 4 + 2], cbuf[i * 4 + 3]);
                }
                loadedReplayColors = true;
            }
        }
    } catch (...) {
    }
    if (!loadedReplayColors) {
        physics::colorFromMass(objs, graphics);
    }

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

    renderer.setProjection(65.0f, 1.0f, 1.0e10f);

    cam.pos = glm::vec3(0.0f, 50.0f, 250.0f);
    cam.front = glm::vec3(0.0f, 0.0f, -1.0f);
    cam.up    = glm::vec3(0.0f, 1.0f,  0.0f);

    renderer.setRenderMode(RenderMode::Sphere);

    Control control(renderer.getWindow(), objs, cam, state);
    control.attach();

    double playbackTime = 0.0;
    double maxTime = (numFrames - 1) * fixedDt; 

    double lastRealTime = glfwGetTime();

    const size_t INVALID_FRAME_IDX = 99999999999ull;
    std::array<std::vector<glm::dvec3>, 4> frames;
    std::array<size_t, 4> frameIds = {INVALID_FRAME_IDX, INVALID_FRAME_IDX, INVALID_FRAME_IDX, INVALID_FRAME_IDX};
    for (auto& f : frames) f.resize(numBodies);
    size_t lastFrameA = INVALID_FRAME_IDX;

    std::vector<glm::dvec3> targetPosition(numBodies);
    std::vector<PosRecord> readScratch(numBodies);
    auto readFrame = [&](size_t frameIdx, std::vector<glm::dvec3>& out) {
        hsize_t count[2] = { 1, static_cast<hsize_t>(numBodies) };
        hsize_t memDims[2] = { 1, static_cast<hsize_t>(numBodies) };
        DataSpace memSpace(2, memDims);

        hsize_t offset[2] = { static_cast<hsize_t>(frameIdx), 0 };
        DataSpace fileSpace = tracksSet.getSpace();
        fileSpace.selectHyperslab(H5S_SELECT_SET, count, offset);
        tracksSet.read(readScratch.data(), GetTrackMemType(), memSpace, fileSpace);
        for (size_t i = 0; i < numBodies; ++i) {
            out[i] = readScratch[i].position;
        }
    };

    double ratio;
    while (!glfwWindowShouldClose(renderer.getWindow()) && state.running) {
        double now = glfwGetTime();
        double realDt = now - lastRealTime;
        lastRealTime = now;
        state.deltaTime = (float)realDt;

        GLFWwindow* win = renderer.getWindow();
        const bool shiftFreeze = (glfwGetKey(win, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
            || (glfwGetKey(win, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS);

        if (!state.pause && !shiftFreeze) {
            playbackTime += realDt * state.timeScale;
        }

        control.updateCameraFromKeys();

        if (playbackTime > maxTime) playbackTime = 0.0;
        if (playbackTime < 0.0) playbackTime = maxTime;
        ratio = playbackTime / fixedDt;
        size_t frameA = (size_t)ratio;
        size_t frameB = frameA + 1;
        if (frameA >= numFrames) frameA = numFrames - 1;
        if (frameB >= numFrames) frameB = numFrames - 1;

        const size_t p0 = (frameA > 0) ? frameA - 1 : frameA;
        const size_t p1 = frameA;
        const size_t p2 = frameB;
        const size_t p3 = std::min(frameB + 1, numFrames - 1);

        const bool movingForward = (lastFrameA != INVALID_FRAME_IDX) && (frameA == lastFrameA + 1);
        if (!movingForward || frameIds[1] != lastFrameA) {
            readFrame(p0, frames[0]); frameIds[0] = p0;
            readFrame(p1, frames[1]); frameIds[1] = p1;
            readFrame(p2, frames[2]); frameIds[2] = p2;
            readFrame(p3, frames[3]); frameIds[3] = p3;
        } else {
            std::swap(frames[0], frames[1]); std::swap(frameIds[0], frameIds[1]);
            std::swap(frames[1], frames[2]); std::swap(frameIds[1], frameIds[2]);
            std::swap(frames[2], frames[3]); std::swap(frameIds[2], frameIds[3]);
            if (frameIds[3] != p3) {
                readFrame(p3, frames[3]);
                frameIds[3] = p3;
            }
        }

        if (frameA < lastFrameA && lastFrameA != INVALID_FRAME_IDX) {
            for (auto& g : graphics) {
                g.resetTrail();
            }
        }

        const double t = ratio - frameA;

        for (size_t i = 0; i < numBodies; ++i) {
            targetPosition[i] = catmullRom(frames[0][i], frames[1][i], frames[2][i], frames[3][i], t);
        }

        for (size_t i = 0; i < numBodies; ++i) {
            objs[i].position = targetPosition[i];
        }

        for (size_t i = 0; i < numBodies; ++i) {
            graphics[i].updateTrail(objs[i].position);
        }

        lastFrameA = frameA;

    renderer.renderFrame(objs, graphics, cam); 

        std::snprintf(title, sizeof(title),
                      "Replay | Time: %.2f / %.2f | Frame: %zu | Sim: x%.2f | Cam: x%.2f%s",
                      playbackTime, maxTime, frameA, state.timeScale, control.getCameraMoveScale(),
                      shiftFreeze ? " | [SHIFT freeze]" : "");
        glfwSetWindowTitle(renderer.getWindow(), title);
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
    
}