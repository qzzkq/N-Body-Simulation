#include "generators.hpp"
#include "physics.hpp"
#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>

static inline float CalcRadiusKm(double m, double rho) {
    const double r_m = std::cbrt((3.0 * m) / (4.0 * 3.14159265358979323846 * rho));
    return static_cast<float>(r_m / 50000.0);
}

class RingScenario : public IScenario {
public:
    std::string getName() const override {
        return "Кольцо из объектов";
    }

    void generate(std::vector<Object>& out, const GenParams& p) override {
        out.clear();
        out.reserve(static_cast<size_t>(p.count) + 1);

        Object center(glm::vec3(0), glm::vec3(0), p.centralMass, 141000.0f, std::nullopt);
        center.Initalizing = false;
        center.radius = CalcRadiusKm(center.mass, center.density);
        out.push_back(center);

        if (p.count <= 0) return;

        std::mt19937 rng(p.seed);
        std::uniform_real_distribution<float> u01(0.0f, 1.0f);
        std::uniform_real_distribution<float> uAngle(0.0f, 6.28318530718f);

        const float ringWidth = std::max(0.0f, p.maxRadius - p.minRadius);
        const float rJitterMax = 0.01f * ringWidth;

        for (int i = 0; i < p.count; ++i) {
            float t = (static_cast<float>(i) + 0.5f) / static_cast<float>(p.count);
            float r2_min = p.minRadius * p.minRadius;
            float r2_max = p.maxRadius * p.maxRadius;
            float r = std::sqrt(r2_min + (r2_max - r2_min) * t);

            r += (u01(rng) * 2.0f - 1.0f) * rJitterMax;
            r = std::clamp(r, p.minRadius, p.maxRadius);

            float a = uAngle(rng);
            float y = (u01(rng) * 2.0f - 1.0f) * (p.spread * 0.5f);

            glm::vec3 pos(r * std::cos(a), y, r * std::sin(a));
            glm::vec3 tdir(-std::sin(a), 0.0f, std::cos(a));

            double dist = std::sqrt(static_cast<double>(r*r) + static_cast<double>(y*y));
            double v_circ_mps = std::sqrt((physics::G * p.centralMass) / (dist * 1000.0));

            float v_kmps = static_cast<float>(v_circ_mps / 1000.0f) * 1.0f;

            Object o(pos, tdir * v_kmps, p.baseMass, 1410.0f, std::nullopt);
            o.Initalizing = false;
            o.radius = CalcRadiusKm(o.mass, o.density);

            out.push_back(std::move(o));
        }

        physics::colorFromMass(out);
    }
};

class ClusterScenario : public IScenario {
public:
    std::string getName() const override {
        return "8 Кластеров";
    }

    void generate(std::vector<Object>& out, const GenParams& p) override {
        out.clear();
        out.reserve(static_cast<size_t>(p.count) + 1);

        Object center(glm::vec3(0), glm::vec3(0), p.centralMass, 141000.0f, std::nullopt);
        center.Initalizing = false;
        center.radius = CalcRadiusKm(center.mass, center.density);
        out.push_back(center);

        if (p.count <= 0) return;

        std::mt19937 rng(p.seed);
        std::uniform_real_distribution<float> u01(0.0f, 1.0f);
        std::uniform_real_distribution<float> uMinus1_1(-1.0f, 1.0f);

        const int numClusters = 8;
        float orbitRadius = (p.minRadius + p.maxRadius) * 0.5f;

        for (int i = 0; i < p.count; ++i) {
            int clusterIdx = i % numClusters;
            float clusterAngle = (2.0f * 3.14159265359f / static_cast<float>(numClusters)) * clusterIdx;

            glm::vec3 clusterCenter(
                orbitRadius * std::cos(clusterAngle),
                0.0f,
                orbitRadius * std::sin(clusterAngle)
            );

            glm::vec3 randomPoint;
            float x, y, z, d2;
            do {
                x = uMinus1_1(rng); y = uMinus1_1(rng); z = uMinus1_1(rng);
                d2 = x*x + y*y + z*z;
            } while (d2 > 1.0f || d2 < 0.0001f);

            float scale = p.spread * std::cbrt(u01(rng)) / std::sqrt(d2);
            glm::vec3 offset(x * scale, y * scale, z * scale);
            glm::vec3 pos = clusterCenter + offset;

            double dist = std::sqrt(static_cast<double>(pos.x*pos.x + pos.z*pos.z));
            float angle = std::atan2(pos.z, pos.x);
            glm::vec3 tdir(-std::sin(angle), 0.0f, std::cos(angle));

            double v_circ_mps = std::sqrt((physics::G * p.centralMass) / (dist * 1000.0));
            float v_kmps = static_cast<float>(v_circ_mps / 1000.0f) * 1.0f;

            Object o(pos, tdir * v_kmps, p.baseMass, 1410.0f, std::nullopt);
            o.Initalizing = false;
            o.radius = CalcRadiusKm(o.mass, o.density);

            out.push_back(std::move(o));
        }
        physics::colorFromMass(out);
    }
};

void ScenarioManager::registerScenario(std::unique_ptr<IScenario> scenario) {
    scenarios_.push_back(std::move(scenario));
}

std::vector<std::string> ScenarioManager::getNames() const {
    std::vector<std::string> names;
    for (const auto& s : scenarios_) {
        names.push_back(s->getName());
    }
    return names;
}

void ScenarioManager::runScenario(size_t index, std::vector<Object>& out, const GenParams& params) {
    if (isValidIndex(index)) {
        std::cout << "Генерируем: " << scenarios_[index]->getName() << "...\n";
        scenarios_[index]->generate(out, params);
    } else {
        std::cerr << "Неправильный индекс генерации!\n";
    }
}

bool ScenarioManager::isValidIndex(size_t index) const {
    return index < scenarios_.size();
}

std::unique_ptr<ScenarioManager> CreateDefaultManager() {
    auto mgr = std::make_unique<ScenarioManager>();
    mgr->registerScenario(std::make_unique<RingScenario>());
    mgr->registerScenario(std::make_unique<ClusterScenario>());
    return mgr;
}