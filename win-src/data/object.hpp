#ifndef OBJECT_HPP
#define OBJECT_HPP

#include <GL/glew.h>
#include <glm/glm.hpp>
#include <optional>
#include <vector>

class Object {
public:
    GLuint VAO, VBO;

    glm::dvec3 position;
    glm::dvec3 velocity;
    glm::vec4 color;

    size_t vertexCount;
    
    bool Initalizing;
    bool Launched;
    bool target;

    double mass;
    double density;
    double radius;
    glm::dvec3 LastPos;

    Object(glm::dvec3 initPosition, glm::dvec3 initVelocity,
           std::optional<double> mass = std::nullopt,
           std::optional<double> density = std::nullopt,
           std::optional<double> radius = std::nullopt);

    void UpdatePos(double deltaTime);
    void UpdateVertices();
    glm::dvec3 GetPos() const;
    void accelerate(double x, double y, double z, double deltaTime);

private:
    std::vector<float> Draw();
};

#endif // OBJECT_HPP