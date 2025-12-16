#pragma once
#include <glm/glm.hpp>
#include <optional>
#include <vector>

class Object {
public:
    // Физические параметры
    glm::dvec3 position;
    glm::dvec3 velocity;
    double mass;
    double density;
    
    // Визуальные параметры
    float radius; 
    glm::vec4 color;

    // Флаги логики
    bool Initalizing = false;
    bool Launched = false;
    bool target = false;

    // Конструктор
    Object(glm::dvec3 initPosition, glm::dvec3 initVelocity,
           std::optional<double> massOpt,
           std::optional<double> densityOpt,
           std::optional<double> radiusOpt);

    // Методы физики
    void UpdatePos(double deltaTime);
    void accelerate(double x, double y, double z, double deltaTime);
    glm::dvec3 GetPos() const;
};