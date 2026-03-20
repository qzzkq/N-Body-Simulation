#pragma once

#include <deque>
#include <glm/glm.hpp>

struct GraphicState {
    glm::vec4 color = glm::vec4(1.0f);
    std::deque<glm::vec3> trail;

    void updateTrail(const glm::dvec3& position, bool force = false) {
        if (force || trail.empty()) {
            trail.push_back(glm::vec3(position));
            return;
        }

        constexpr float distanceThreshold = 0.2f;
        if (glm::distance(glm::vec3(position), trail.back()) > distanceThreshold) {
            trail.push_back(glm::vec3(position));
        }

        if (trail.size() > MAX_TRAIL_LENGTH) {
            trail.pop_front();
        }
    }

    void resetTrail() {
        trail.clear();
    }

    static constexpr std::size_t MAX_TRAIL_LENGTH = 150;
};
