#include "renderer.hpp"
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>
#include <cmath>
#include <algorithm>

static const char* VERTEX_SHADER_SOURCE = R"glsl(
#version 330 core
layout(location=0) in vec3 aPos;
layout(location=2) in vec4 iModel0;
layout(location=3) in vec4 iModel1;
layout(location=4) in vec4 iModel2;
layout(location=5) in vec4 iModel3;
layout(location=6) in vec4 iColor;
uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform vec4 objectColor;
uniform int useInstancing;
uniform float farPlane;
out vec4 vColor;
void main() {
    mat4 instanceMatrix = mat4(iModel0, iModel1, iModel2, iModel3);
    mat4 M = (useInstancing == 1) ? instanceMatrix : model;
    vColor = (useInstancing == 1) ? iColor : objectColor;
    gl_Position = projection * view * M * vec4(aPos, 1.0);
    float C = 1.0;
    gl_Position.z = (2.0 * log(C * gl_Position.w + 1.0) / log(C * farPlane + 1.0) - 1.0) * gl_Position.w;
}
)glsl";

static const char* FRAGMENT_SHADER_SOURCE = R"glsl(
#version 330 core
out vec4 FragColor;
in vec4 vColor;
void main() {
    FragColor = vColor;
}
)glsl";

static const char* TRAIL_VS = R"glsl(
#version 330 core
layout(location=0) in vec3 aPos;
uniform mat4 view;
uniform mat4 projection;
void main() {
    gl_Position = projection * view * vec4(aPos, 1.0);
}
)glsl";

static const char* TRAIL_FS = R"glsl(
#version 330 core
out vec4 FragColor;
uniform vec4 color;
void main() {
    FragColor = color;
}
)glsl";

static const char* POINT_VS = R"glsl(
#version 330 core
layout(location=0) in vec3 aPos;
layout(location=1) in vec4 aColor;
uniform mat4 view;
uniform mat4 projection;
uniform float farPlane;
out vec4 vColor;
void main() {
    vColor = aColor;
    gl_Position = projection * view * vec4(aPos, 1.0);
    float C = 1.0;
    gl_Position.z = (2.0 * log(C * gl_Position.w + 1.0) / log(C * farPlane + 1.0) - 1.0) * gl_Position.w;
    float dist = length(aPos);
    const float refDist = 2000.0;
    const float pxRef = 4.5;
    float sizePx = pxRef * refDist / max(dist, refDist * 0.02);
    gl_PointSize = clamp(sizePx, 0.5, 48.0);
}
)glsl";

static const char* POINT_FS = R"glsl(
#version 330 core
in vec4 vColor;
out vec4 FragColor;
void main() {
    FragColor = vColor;
}
)glsl";

// перевод сферических координат в декартовы 
static glm::vec3 sphericalToCartesian(float r, float theta, float phi) {
    return glm::vec3(
        r * sin(theta) * cos(phi),
        r * cos(theta),
        r * sin(theta) * sin(phi)
    );
}

//==PUBLIC BLOCK==

Renderer::Renderer() {
}

Renderer::~Renderer() {
    if (successInit_) {
        glDeleteVertexArrays(1, &sphereVAO_); glDeleteBuffers(1, &sphereVBO_);
        glDeleteVertexArrays(1, &cubeVAO_);   glDeleteBuffers(1, &cubeVBO_);
        glDeleteVertexArrays(1, &pointVAO_);
        glDeleteBuffers(1, &pointVBO_);
        glDeleteBuffers(1, &pointColorVBO_);
        glDeleteBuffers(1, &instanceModelVBO_);
        glDeleteBuffers(1, &instanceColorVBO_);
        glDeleteProgram(program_);
        glDeleteProgram(pointProgram_);
    }
}

bool Renderer::init(int width,  int height, const char* title, bool fullscreen, bool maximized) {
    if (!initWindow(width, height, title, fullscreen, maximized)) {
        successInit_ = false; 
        return false; 
    }

    initProgram();

    initCubeGeometry();
    initPointGeometry();
    initSphereGeometry();
    initTrailVAO_VBO(); 
    successInit_ = true;  
    return true; 
}

void Renderer::setProjection(float fov_deg,
                             float znear,
                             float zfar)
{
    int width, height;
    glfwGetFramebufferSize(window_, &width, &height);
    float aspect = (float) width / height; 
    farPlane_ = zfar;
    glUseProgram(program_);
    projectionMatrix_ = glm::perspective(glm::radians(fov_deg), aspect, znear, zfar);
    glUniformMatrix4fv(uProj_, 1, GL_FALSE, glm::value_ptr(projectionMatrix_));
    glUniform1f(glGetUniformLocation(program_, "farPlane"), farPlane_);
    glUseProgram(pointProgram_);
    glUniformMatrix4fv(glGetUniformLocation(pointProgram_, "projection"), 1, GL_FALSE, glm::value_ptr(projectionMatrix_));
    glUniform1f(glGetUniformLocation(pointProgram_, "farPlane"), farPlane_);
    glUseProgram(trailProgram_);
    glUniformMatrix4fv(glGetUniformLocation(trailProgram_, "projection"), 1, GL_FALSE, glm::value_ptr(projectionMatrix_));
}

void Renderer::updateView(const Camera& camera)
{
    // Relative-to-eye рендеринг: translation already baked into model positions.
    // View matrix contains orientation only.
    viewMatrix_ = glm::lookAt(glm::vec3(0.0f), camera.front, camera.up);
    glUseProgram(program_);
    glUniformMatrix4fv(uView_, 1, GL_FALSE, glm::value_ptr(viewMatrix_));
    glUseProgram(pointProgram_);
    glUniformMatrix4fv(glGetUniformLocation(pointProgram_, "view"), 1, GL_FALSE, glm::value_ptr(viewMatrix_));
    glUseProgram(trailProgram_);
    glUniformMatrix4fv(glGetUniformLocation(trailProgram_, "view"), 1, GL_FALSE, glm::value_ptr(viewMatrix_));
}

void Renderer::drawObjects(const std::vector<Object>& objs, const std::vector<GraphicState>& graphics, const Camera& cam) {
    glUseProgram(program_);

    if (graphics.size() != objs.size()) return;

    static std::vector<glm::mat4> modelMatrices;
    static std::vector<glm::vec4> colors;
    static std::vector<glm::vec3> pointPositions;
    static std::vector<glm::vec4> pointColors;

    modelMatrices.clear(); colors.clear();
    pointPositions.clear(); pointColors.clear();

    const float lodSphereMaxDist = 15000.0f; 

    for (std::size_t i = 0; i < objs.size(); ++i) {
        glm::dvec3 relativePos = objs[i].position - glm::dvec3(cam.pos);
        float distance = static_cast<float>(glm::length(relativePos));

        // Если режим Точки ИЛИ объект далеко — рисуем точкой
        if (mode_ == RenderMode::Points || distance >= lodSphereMaxDist) {
            pointPositions.push_back(glm::vec3(relativePos));
            pointColors.push_back(graphics[i].color);
        } else {
            // Иначе (Сферы или Кубы вблизи)
            const float farFalloff = 1.0f / (1.0f + distance / 30000.0f);
            float visualScale = std::max(static_cast<float>(objs[i].radius), distance * 0.002f * farFalloff);

            glm::mat4 M(1.0f);
            M = glm::translate(M, glm::vec3(relativePos));
            M = glm::scale(M, glm::vec3(visualScale));

            modelMatrices.push_back(M);
            colors.push_back(graphics[i].color);
        }
    }

    // Отрисовка 3D-геометрии (Кубы или Сферы)
    if (!modelMatrices.empty()) {
        glBindBuffer(GL_ARRAY_BUFFER, instanceModelVBO_);
        glBufferData(GL_ARRAY_BUFFER, modelMatrices.size() * sizeof(glm::mat4), modelMatrices.data(), GL_DYNAMIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, instanceColorVBO_);
        glBufferData(GL_ARRAY_BUFFER, colors.size() * sizeof(glm::vec4), colors.data(), GL_DYNAMIC_DRAW);

        glUniform1i(uUseInstancing_, 1);
        
        if (mode_ == RenderMode::Cubes) {
            glBindVertexArray(cubeVAO_);
            glDrawArraysInstanced(GL_TRIANGLES, 0, 36, static_cast<GLsizei>(modelMatrices.size()));
        } else {
            glBindVertexArray(sphereVAO_);
            glDrawArraysInstanced(GL_TRIANGLES, 0, static_cast<GLsizei>(sphereVertexCount_), static_cast<GLsizei>(modelMatrices.size()));
        }
        glBindVertexArray(0);
    }

    //Отрисовка точек
    if (!pointPositions.empty()) {
        glUseProgram(pointProgram_);
        glBindBuffer(GL_ARRAY_BUFFER, pointVBO_);
        glBufferData(GL_ARRAY_BUFFER, pointPositions.size() * sizeof(glm::vec3), pointPositions.data(), GL_DYNAMIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, pointColorVBO_);
        glBufferData(GL_ARRAY_BUFFER, pointColors.size() * sizeof(glm::vec4), pointColors.data(), GL_DYNAMIC_DRAW);

        glBindVertexArray(pointVAO_);
        glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(pointPositions.size()));
        glBindVertexArray(0);
    }

    glUseProgram(program_);
}

void Renderer::resizeWindow(int w, int h) {
    glfwSetWindowSize(window_, w, h);
    int frW, frH;
    glfwGetFramebufferSize(window_, &frW, &frH);
    glViewport(0, 0, frW, frH);
}

void Renderer::renderFrame(const std::vector<Object>& objs, const std::vector<GraphicState>& graphics, const Camera& cam) {
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    updateView(cam);
    
    drawTrails(objs, graphics, cam);
    drawObjects(objs, graphics, cam);
    
    glfwSwapBuffers(window_);
}

void Renderer::drawTrails(const std::vector<Object>& objs, const std::vector<GraphicState>& graphics, const Camera& cam) {
    if (graphics.size() != objs.size()) return;

    glUseProgram(trailProgram_);

    glBindVertexArray(trailVAO_);
    glBindBuffer(GL_ARRAY_BUFFER, trailVBO_);

    // std::deque не гарантирует непрерывное расположение в памяти, поэтому
    // перед отправкой в GPU копируем в небольшой временный буфер.
    static thread_local std::vector<glm::vec3> trailScratch;
    for (std::size_t i = 0; i < objs.size(); ++i) {
        const auto& gfx = graphics[i];
        if (gfx.trail.size() < 2) continue;

        trailScratch.resize(gfx.trail.size());
        for (std::size_t j = 0; j < gfx.trail.size(); ++j) {
            glm::dvec3 relativeTrailPos = glm::dvec3(gfx.trail[j]) - glm::dvec3(cam.pos);
            trailScratch[j] = glm::vec3(relativeTrailPos);
        }
        glBufferData(GL_ARRAY_BUFFER, gfx.trail.size() * sizeof(glm::vec3), trailScratch.data(), GL_DYNAMIC_DRAW);
        
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
        glEnableVertexAttribArray(0);

        glm::vec4 trailColor = gfx.color * 0.7f;
        trailColor.a = 0.5f; 
        glUniform4fv(glGetUniformLocation(trailProgram_, "color"), 1, glm::value_ptr(trailColor));

        glDrawArrays(GL_LINE_STRIP, 0, static_cast<GLsizei>(gfx.trail.size()));
    }
    glBindVertexArray(0);
}

//==PRIVATE BLOCK==

// Private-методы для компиляции шефдерной программы 

GLuint Renderer::compileProgram(const char* vs, const char* fs) {
    auto compile = [](GLenum type, const char* src){
        GLuint s = glCreateShader(type);
        glShaderSource(s, 1, &src, nullptr);
        glCompileShader(s);
        GLint ok=0; glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
        if(!ok){
            char log[1024];
            glGetShaderInfoLog(s,1024,nullptr,log);
            std::cerr << "Shader compile error: " << log << std::endl;
        }
        return s;
    };
    GLuint v = compile(GL_VERTEX_SHADER,   vs);
    GLuint f = compile(GL_FRAGMENT_SHADER, fs);
    GLuint p = glCreateProgram();
    glAttachShader(p, v);
    glAttachShader(p, f);
    glLinkProgram(p);
    glDeleteShader(v);
    glDeleteShader(f);
    return p;
}

// Private-методы для инициализации отрисовочных параметров 

void Renderer::initPointGeometry() {
    glGenVertexArrays(1, &pointVAO_);
    glGenBuffers(1, &pointVBO_);
    glGenBuffers(1, &pointColorVBO_);

    glBindVertexArray(pointVAO_);

    // Позиции точек мира.
    glBindBuffer(GL_ARRAY_BUFFER, pointVBO_);
    glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // Цвета точек.
    glBindBuffer(GL_ARRAY_BUFFER, pointColorVBO_);
    glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1);

    glBindVertexArray(0);
}

void Renderer::initCubeGeometry() {
    // Стандартный куб -1..1
    float cubeVerts[] = {
        // back
        -1.f,-1.f,-1.f,  1.f, 1.f,-1.f,  1.f,-1.f,-1.f,
        -1.f,-1.f,-1.f, -1.f, 1.f,-1.f,  1.f, 1.f,-1.f,
        // front
        -1.f,-1.f, 1.f,  1.f,-1.f, 1.f,  1.f, 1.f, 1.f,
        -1.f,-1.f, 1.f,  1.f, 1.f, 1.f, -1.f, 1.f, 1.f,
        // left
        -1.f,-1.f,-1.f, -1.f,-1.f, 1.f, -1.f, 1.f, 1.f,
        -1.f,-1.f,-1.f, -1.f, 1.f, 1.f, -1.f, 1.f,-1.f,
        // right
         1.f,-1.f,-1.f,  1.f, 1.f, 1.f,  1.f,-1.f, 1.f,
         1.f,-1.f,-1.f,  1.f, 1.f,-1.f,  1.f, 1.f, 1.f,
        // bottom
        -1.f,-1.f,-1.f,  1.f,-1.f, 1.f,  1.f,-1.f,-1.f,
        -1.f,-1.f,-1.f, -1.f,-1.f, 1.f,  1.f,-1.f, 1.f,
        // top
        -1.f, 1.f,-1.f,  1.f, 1.f,-1.f,  1.f, 1.f, 1.f,
        -1.f, 1.f,-1.f,  1.f, 1.f, 1.f, -1.f, 1.f, 1.f
    };
    glGenVertexArrays(1, &cubeVAO_);
    glGenBuffers(1, &cubeVBO_);
    glBindVertexArray(cubeVAO_);
    glBindBuffer(GL_ARRAY_BUFFER, cubeVBO_);
    glBufferData(GL_ARRAY_BUFFER, sizeof(cubeVerts), cubeVerts, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    setupInstanceAttribs(cubeVAO_);
    glBindVertexArray(0);
}

void Renderer::initSphereGeometry() {
    std::vector<float> vertices;
    const int stacks = 12;  // Детализация
    const int sectors = 12;
    const float radius = 1.0f; // Единичная сфера

    for (int i = 0; i < stacks; ++i) {
        float theta1 = (float(i) / stacks) * glm::pi<float>();
        float theta2 = (float(i + 1) / stacks) * glm::pi<float>();

        for (int j = 0; j < sectors; ++j) {
            float phi1 = (float(j) / sectors) * 2.0f * glm::pi<float>();
            float phi2 = (float(j + 1) / sectors) * 2.0f * glm::pi<float>();

            glm::vec3 v1 = sphericalToCartesian(radius, theta1, phi1);
            glm::vec3 v2 = sphericalToCartesian(radius, theta1, phi2);
            glm::vec3 v3 = sphericalToCartesian(radius, theta2, phi1);
            glm::vec3 v4 = sphericalToCartesian(radius, theta2, phi2);

            // Два треугольника на сектор
            auto addV = [&](const glm::vec3& v) {
                vertices.push_back(v.x); vertices.push_back(v.y); vertices.push_back(v.z);
            };

            addV(v1); addV(v2); addV(v3);
            addV(v2); addV(v4); addV(v3);
        }
    }
    
    sphereVertexCount_ = vertices.size() / 3;

    glGenVertexArrays(1, &sphereVAO_);
    glGenBuffers(1, &sphereVBO_);
    
    glBindVertexArray(sphereVAO_);
    glBindBuffer(GL_ARRAY_BUFFER, sphereVBO_);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);
    
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    setupInstanceAttribs(sphereVAO_);
    glBindVertexArray(0);
}

void Renderer::initTrailVAO_VBO() {
    glGenVertexArrays(1, &trailVAO_);
    glGenBuffers(1, &trailVBO_);
}

bool Renderer::initWindow(int width, int height, const char* title, bool fullscreen, bool maximized) {
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW\n";
        return false;
    }
    GLFWmonitor* primaryMonitor = glfwGetPrimaryMonitor();
    const GLFWvidmode* mode = glfwGetVideoMode(primaryMonitor);

    if (mode == nullptr) {
        std::cerr << "Failed to get video mode\n";
        glfwTerminate();
        return false;
    }

    GLFWmonitor* monitorForWindow = nullptr; 

    if (fullscreen) {
        monitorForWindow = primaryMonitor;
        width = mode->width;   
        height = mode->height;
        glfwWindowHint(GLFW_REFRESH_RATE, mode->refreshRate);
    } 
    else if (maximized) {
        glfwWindowHint(GLFW_MAXIMIZED, GLFW_TRUE);
        
        width = mode->width;
        height = mode->height;
    }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow(width, height, title, monitorForWindow, nullptr);
    
    if (!window) {
        std::cerr << "Failed to create GLFW window.\n";
        glfwTerminate();
        return false;
    }
    glfwMakeContextCurrent(window);

    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW.\n";
        glfwTerminate();
        return false;
    }

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_PROGRAM_POINT_SIZE);

    int bufferWidth, bufferHeight;
    glfwGetFramebufferSize(window, &bufferWidth, &bufferHeight);
    glViewport(0, 0, bufferWidth, bufferHeight);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    window_ = window; 
    return true;
}

void Renderer::initProgram() {
    program_ = compileProgram(VERTEX_SHADER_SOURCE, FRAGMENT_SHADER_SOURCE);
    trailProgram_ = compileProgram(TRAIL_VS, TRAIL_FS);
    pointProgram_ = compileProgram(POINT_VS, POINT_FS);
    glUseProgram(program_);
    
    uModel_ = glGetUniformLocation(program_, "model");
    uView_  = glGetUniformLocation(program_, "view");
    uProj_  = glGetUniformLocation(program_, "projection");
    uColor_ = glGetUniformLocation(program_, "objectColor");
    uUseInstancing_ = glGetUniformLocation(program_, "useInstancing");

    glGenBuffers(1, &instanceModelVBO_);
    glGenBuffers(1, &instanceColorVBO_);
    glUniform1f(glGetUniformLocation(program_, "farPlane"), farPlane_);

    glUseProgram(pointProgram_);
    glUniform1f(glGetUniformLocation(pointProgram_, "farPlane"), farPlane_);
}

void Renderer::setupInstanceAttribs(GLuint vao) {
    glBindVertexArray(vao);

    glBindBuffer(GL_ARRAY_BUFFER, instanceModelVBO_);
    constexpr GLsizei matStride = static_cast<GLsizei>(sizeof(glm::mat4));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, matStride, (void*)0);
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, matStride, (void*)(sizeof(glm::vec4)));
    glEnableVertexAttribArray(4);
    glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, matStride, (void*)(2 * sizeof(glm::vec4)));
    glEnableVertexAttribArray(5);
    glVertexAttribPointer(5, 4, GL_FLOAT, GL_FALSE, matStride, (void*)(3 * sizeof(glm::vec4)));
    glVertexAttribDivisor(2, 1);
    glVertexAttribDivisor(3, 1);
    glVertexAttribDivisor(4, 1);
    glVertexAttribDivisor(5, 1);

    glBindBuffer(GL_ARRAY_BUFFER, instanceColorVBO_);
    glEnableVertexAttribArray(6);
    glVertexAttribPointer(6, 4, GL_FLOAT, GL_FALSE, sizeof(glm::vec4), (void*)0);
    glVertexAttribDivisor(6, 1);

    glBindVertexArray(0);
}