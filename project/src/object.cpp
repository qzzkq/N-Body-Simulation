#include "object.hpp"
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <cmath>

static glm::dvec3 sphericalToCartesian(double r, double theta, double phi) {
    double x = r * sin(theta) * cos(phi);
    double y = r * cos(theta);
    double z = r * sin(theta) * sin(phi);
    return glm::dvec3(x, y, z);
}

Object::Object(glm::dvec3 initPosition, glm::dvec3 initVelocity, // конструктор - установка параметров объекта 
               std::optional<double> massOpt,
               std::optional<double> densityOpt,
               std::optional<double> radiusOpt)
{
    this->position = initPosition;
    this->velocity = initVelocity;
    this->color = glm::vec4(1.0f, 0.0f, 0.0f, 1.0f);

    // Вычисление отсутствующего параметра при необходимости
    if (massOpt.has_value() && densityOpt.has_value() && radiusOpt.has_value()) {
        this->mass = *massOpt;
        this->density = *densityOpt;
        this->radius = *radiusOpt;
    } else if (massOpt.has_value() && densityOpt.has_value()) {
        this->mass    = *massOpt;
        this->density = *densityOpt;
        this->radius  = std::pow((3.0 * mass / density) / (4.0 * 3.14159265359), 1.0/3.0) / 100000.0;
    } else if (massOpt.has_value() && radiusOpt.has_value()) {
        double volume = (4.0 / 3.0) * 3.14159265359 * std::pow(*radiusOpt, 3);
        double dens   = *massOpt / volume;
        this->density = dens;
        this->mass    = *massOpt;
        this->radius  = *radiusOpt;
    } else if (radiusOpt.has_value() && densityOpt.has_value()) {
        double volume = (4.0 / 3.0) * 3.14159265359 * std::pow(*radiusOpt, 3);
        double m = *densityOpt * volume;
        this->mass    = m;
        this->density = *densityOpt;
        this->radius  = *radiusOpt;
    }
}

void Object::UpdatePos(double deltaTime) {
    this->position += velocity * deltaTime;
    this->radius = std::pow((3.0 * mass / density) / (4.0 * 3.14159265359), 1.0/3.0) / 100000.0;
}

// геттеры и сеттеры 
glm::dvec3 Object::GetPos() const {
    return this->position;
}
void Object::SetPos(glm::dvec3 newPos) {
    this->position = newPos; 
}

glm::dvec3 Object::GetVel() const {
    return this->velocity; 
}
void Object::SetVel(glm::dvec3 newVel) {
    this->velocity = newVel; 
}

glm::vec4 Object::GetColor() const {
    return this->color; 
}
void Object::SetColor(glm::vec4 newColor) {
    this->color = newColor; 
}

double Object::GetMass() const {
    return this->mass; 
}
void Object::SetMass(double newMass) {
    this->mass = newMass; 
} 

double Object::GetDensity() const {
    return this->density; 
}
void Object::SetDensity(double newDensity) {
    this->density = newDensity; 
}

double Object::GetRadius() const {
    return this->radius; 
}
void Object::SetRadius(double newRadius) {
    this->radius = newRadius;
}


void Object::accelerate(double x, double y, double z, double deltaTime) {
    this->velocity.x += x * deltaTime;
    this->velocity.y += y * deltaTime;
    this->velocity.z += z * deltaTime;
}
