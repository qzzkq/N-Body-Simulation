#include "object.hpp"
#include <cmath>
#include <algorithm>

Object::Object(glm::dvec3 initPosition, glm::dvec3 initVelocity,
               std::optional<double> massOpt,
               std::optional<double> densityOpt,
               std::optional<double> radiusOpt)
{
    position = initPosition;
    velocity = initVelocity;
    Initalizing = false;
    Launched = false;
    color = glm::vec4(1.0f); // Белый по умолчанию

    // Логика вычисления недостающих параметров
    if (massOpt.has_value() && densityOpt.has_value() && radiusOpt.has_value()) {
        mass = *massOpt; density = *densityOpt; radius = static_cast<float>(*radiusOpt);
    } 
    else if (massOpt.has_value() && densityOpt.has_value()) {
        mass = *massOpt; 
        density = *densityOpt;
        double r = std::pow((3.0 * mass / density) / (4.0 * 3.14159265359), 1.0/3.0);
        radius = static_cast<float>(r / 50000.0); 
    } 
    else if (massOpt.has_value() && radiusOpt.has_value()) {
        mass = *massOpt; 
        radius = static_cast<float>(*radiusOpt);
        double vol = (4.0 / 3.0) * 3.14159265359 * std::pow(radius, 3);
        density = mass / vol;
    }
    
    else {
        mass = 1.0; density = 1.0; radius = 1.0f;
    }
}

void Object::UpdatePos(double deltaTime) {
    position += velocity * deltaTime;
}

void Object::accelerate(double x, double y, double z, double deltaTime) {
    velocity.x += x * deltaTime;
    velocity.y += y * deltaTime;
    velocity.z += z * deltaTime;
}

glm::dvec3 Object::GetPos() const {
    return position;
}
