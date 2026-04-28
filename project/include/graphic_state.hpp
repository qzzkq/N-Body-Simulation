#pragma once

#include <array>
#include <glm/glm.hpp>

static constexpr std::size_t MAX_TRAIL_LENGTH = 150;

struct GraphicState {
    glm::vec4 color = glm::vec4(1.0f);

    std::array<glm::vec3, MAX_TRAIL_LENGTH> trail{};
    int trailHead = 0;  // ring buffer write pointer
    int trailSize = 0;  // [0, MAX_TRAIL_LENGTH]

    // Логический индекс: 0 = самая старая точка, trailSize-1 = самая новая
    const glm::vec3& trailAt(int i) const {
        int start = (trailHead - trailSize + static_cast<int>(MAX_TRAIL_LENGTH)) % static_cast<int>(MAX_TRAIL_LENGTH);
        return trail[static_cast<std::size_t>((start + i) % static_cast<int>(MAX_TRAIL_LENGTH))];
    }

    void updateTrail(const glm::dvec3& position, bool force = false) {
        const glm::vec3 pos32(position);

        if (!force && trailSize > 0) {
            constexpr float distanceThreshold = 0.2f;
            const glm::vec3& last = trail[static_cast<std::size_t>(
                (trailHead - 1 + static_cast<int>(MAX_TRAIL_LENGTH)) % static_cast<int>(MAX_TRAIL_LENGTH))];
            if (glm::distance(pos32, last) <= distanceThreshold) return;
        }

        trail[static_cast<std::size_t>(trailHead)] = pos32;
        trailHead = (trailHead + 1) % static_cast<int>(MAX_TRAIL_LENGTH);
        if (trailSize < static_cast<int>(MAX_TRAIL_LENGTH)) ++trailSize;
    }

    void resetTrail() {
        trailHead = 0;
        trailSize = 0;
    }
};
