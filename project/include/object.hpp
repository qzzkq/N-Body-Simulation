#ifndef OBJECT_HPP
#define OBJECT_HPP

#include <optional>
#include <glm/glm.hpp>
#include <vector>

class Object {
public:

    Object(glm::dvec3 initPosition, glm::dvec3 initVelocity,
           std::optional<double> mass = std::nullopt,
           std::optional<double> density = std::nullopt,
           std::optional<double> radius = std::nullopt);

    void UpdatePos(double deltaTime);
    void accelerate(double x, double y, double z, double deltaTime);

    // геттеры и сеттеры 
    glm::dvec3 GetPos() const;
    void SetPos(glm::dvec3 newPos); 

    glm::dvec3 GetVel() const;
    void SetVel(glm::dvec3 newVel); 

    glm::vec4 GetColor() const;
    void SetColor(glm::vec4 newColor); 

    double GetMass() const;
    void SetMass(double newMass); 

    double GetDensity() const;
    void SetDensity(double newDensity); 

    double GetRadius() const;
    void SetRadius(double newRadius); 

private:
    glm::dvec3 position;
    glm::dvec3 velocity;
    glm::vec4 color;
    double mass;
    double density;
    double radius;
};

#endif // OBJECT_HPP