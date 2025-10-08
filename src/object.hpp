#ifndef OBJECT_HPP
#define OBJECT_HPP

#include <GL/glew.h>
#include <glm/glm.hpp>
#include <optional>
#include <vector>

class Object {
public:
    GLuint VAO, VBO;

    glm::vec3 position;
    glm::vec3 velocity;
    glm::vec4 color;

    size_t vertexCount;
    
    bool Initalizing;
    bool Launched;
    bool target;

    float mass;
    float density;
    float radius;
    glm::vec3 LastPos;

    Object(glm::vec3 initPosition, glm::vec3 initVelocity,
           std::optional<float> mass = std::nullopt,
           std::optional<float> density = std::nullopt,
           std::optional<float> radius = std::nullopt);

    void UpdatePos(float deltaTime);
    void UpdateVertices();
    glm::vec3 GetPos() const;
    void accelerate(float x, float y, float z, float deltaTime);

private:
    std::vector<float> Draw();
};

#endif // OBJECT_HPP