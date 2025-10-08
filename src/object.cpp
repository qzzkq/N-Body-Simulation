#include "object.hpp"
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/gtc/constants.hpp>
#include <cmath>

static glm::vec3 sphericalToCartesian(float r, float theta, float phi) {
    float x = r * sin(theta) * cos(phi);
    float y = r * cos(theta);
    float z = r * sin(theta) * sin(phi);
    return glm::vec3(x, y, z);
}

Object::Object(glm::vec3 initPosition, glm::vec3 initVelocity,
               std::optional<float> massOpt,
               std::optional<float> densityOpt,
               std::optional<float> radiusOpt)
{
    position = initPosition;
    velocity = initVelocity;
    Initalizing = false;
    Launched = false;
    target = false;
    color = glm::vec4(1.0f, 0.0f, 0.0f, 1.0f);
    LastPos = position;

    // Вычисление отсутствующего параметра при необходимости
    if (massOpt.has_value() && densityOpt.has_value() && radiusOpt.has_value()) {
        mass = *massOpt;
        density = *densityOpt;
        radius = *radiusOpt;
    } else if (massOpt.has_value() && densityOpt.has_value()) {
        mass    = *massOpt;
        density = *densityOpt;
        radius  = std::pow((3 * mass / density) / (4 * 3.14159265359f), 1.0f/3.0f) / 100000.0f;
    } else if (massOpt.has_value() && radiusOpt.has_value()) {
        float volume = (4.0f / 3.0f) * 3.14159265359f * std::pow(*radiusOpt, 3);
        float dens   = *massOpt / volume;
        density = dens;
        mass    = *massOpt;
        radius  = *radiusOpt;
    } else if (radiusOpt.has_value() && densityOpt.has_value()) {
        float volume = (4.0f / 3.0f) * 3.14159265359f * std::pow(*radiusOpt, 3);
        float m = *densityOpt * volume;
        mass    = m;
        density = *densityOpt;
        radius  = *radiusOpt;
    }

    std::vector<float> vertices = Draw();
    vertexCount = vertices.size();

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);
}

// Генерация вершин для сферы 
std::vector<float> Object::Draw() {
    std::vector<float> vertices;
    int stacks = 10;
    int sectors = 10;
    for (int i = 0; i <= stacks; ++i) {
        float theta1 = (i / (float)stacks) * glm::pi<float>();
        float theta2 = ((i+1) / (float)stacks) * glm::pi<float>();
        for (int j = 0; j < sectors; ++j) {
            float phi1 = (j / (float)sectors) * 2 * glm::pi<float>();
            float phi2 = ((j+1) / (float)sectors) * 2 * glm::pi<float>();

            glm::vec3 v1 = sphericalToCartesian(radius, theta1, phi1);
            glm::vec3 v2 = sphericalToCartesian(radius, theta1, phi2);
            glm::vec3 v3 = sphericalToCartesian(radius, theta2, phi1);
            glm::vec3 v4 = sphericalToCartesian(radius, theta2, phi2);

            // Первый треугольник 
            vertices.insert(vertices.end(), {v1.x, v1.y, v1.z});
            vertices.insert(vertices.end(), {v2.x, v2.y, v2.z});
            vertices.insert(vertices.end(), {v3.x, v3.y, v3.z});
            // Второй треугольник 
            vertices.insert(vertices.end(), {v2.x, v2.y, v2.z});
            vertices.insert(vertices.end(), {v4.x, v4.y, v4.z});
            vertices.insert(vertices.end(), {v3.x, v3.y, v3.z});
        }
    }
    return vertices;
}

void Object::UpdatePos(float deltaTime) {
    position += velocity * deltaTime;
    radius = std::pow((3 * mass / density) / (4 * 3.14159265359f), 1.0f/3.0f) / 100000.0f;
}

void Object::UpdateVertices() {
    std::vector<float> vertices = Draw();
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);
}

glm::vec3 Object::GetPos() const {
    return position;
}

void Object::accelerate(float x, float y, float z, float deltaTime) {
    velocity.x += x * deltaTime;
    velocity.y += y * deltaTime;
    velocity.z += z * deltaTime;
}
