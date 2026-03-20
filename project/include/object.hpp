#pragma once
#include <glm/glm.hpp>
#include <optional>
#include <string>

class Object {
public:
    glm::dvec3 position;
    glm::dvec3 velocity;
    glm::dvec3 acceleration = glm::dvec3(0.0);
    double mass;
    double density;
    
    float radius; 
    std::string name;

    bool Initalizing = false;
    bool Launched = false;
    bool target = false;

    Object(glm::dvec3 initPosition, glm::dvec3 initVelocity,
           std::optional<double> massOpt,
           std::optional<double> densityOpt,
           std::optional<double> radiusOpt);

    void UpdatePos(double deltaTime);
    void accelerate(double x, double y, double z, double deltaTime);
    glm::dvec3 GetPos() const;
};