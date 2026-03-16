#pragma once
#include <glm/glm.hpp>
#include <optional>
#include <vector>
#include <deque> 
#include <string>

class Object {
public:
    glm::dvec3 position;
    glm::dvec3 velocity;
    double mass;
    double density;
    
    float radius; 
    glm::vec4 color;
    std::string name;
    std::vector<glm::vec3> trail;

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

    void updateTrail(bool force = false) {
        if (force || trail.empty()) { 
            trail.push_back(position);
            return;
        }

        float distanceThreshold = 1.0f;
        if (glm::distance(glm::vec3(position), trail.back()) > distanceThreshold) {
            trail.push_back(position);
        }

        if (trail.size() > MAX_TRAIL_LENGTH) {
            trail.erase(trail.begin());
        }
    }

    void resetTrail() {
        trail.clear();
    }

private: 
    static const size_t MAX_TRAIL_LENGTH = 150; 
};