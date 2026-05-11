#include "generators.hpp"
#include "physics.hpp"
#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>

static inline float CalcRadiusAU(double m) {
    const double clampedMass = std::max(m, 1e-12);
    return static_cast<float>(0.02 * std::cbrt(clampedMass));
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
        center.radius = CalcRadiusAU(center.mass);
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
            double v_circ_au_per_year = std::sqrt((physics::G * p.centralMass) / std::max(dist, 1e-9));

            Object o(pos, tdir * static_cast<float>(v_circ_au_per_year), p.baseMass, 1410.0f, std::nullopt);
            o.Initalizing = false;
            o.radius = CalcRadiusAU(o.mass);

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
        center.radius = CalcRadiusAU(center.mass);
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

            double v_circ_au_per_year = std::sqrt((physics::G * p.centralMass) / std::max(dist, 1e-9));

            Object o(pos, tdir * static_cast<float>(v_circ_au_per_year), p.baseMass, 1410.0f, std::nullopt);
            o.Initalizing = false;
            o.radius = CalcRadiusAU(o.mass);

            out.push_back(std::move(o));
        }
        physics::colorFromMass(out);
    }
};

// Шаровое скопление по модели Plummer (заготовка для M13).
// Total mass и масштабный радиус берутся из GenParams: centralMass = M_total [M☉],
// maxRadius = Plummer scale a [AU]; при нулях используются значения для M13.
class M13PlummerScenario : public IScenario {
public:
    std::string getName() const override {
        return "M13 (Plummer-скопление)";
    }

    void generate(std::vector<Object>& out, const GenParams& p) override {
        out.clear();
        if (p.count <= 0) return;
        out.reserve(static_cast<size_t>(p.count));

        const double M_total_solar = (p.centralMass > 0.0) ? p.centralMass : 6.0e5;
        const double a_au          = (p.maxRadius   > 0.0) ? static_cast<double>(p.maxRadius) : 8000.0;
        const double m_star        = M_total_solar / static_cast<double>(p.count);
        const double GM            = physics::G * M_total_solar;
        const double a2            = a_au * a_au;

        std::mt19937 rng(p.seed);
        std::uniform_real_distribution<double> u01(0.0, 1.0);
        std::uniform_real_distribution<double> uMinus1_1(-1.0, 1.0);
        std::uniform_real_distribution<double> uAngle(0.0, 2.0 * 3.14159265358979323846);

        for (int i = 0; i < p.count; ++i) {
            // Plummer CDF: r = a / sqrt(u^{-2/3} - 1).
            double u = u01(rng);
            if (u < 1e-12) u = 1e-12;
            const double r = a_au / std::sqrt(std::pow(u, -2.0 / 3.0) - 1.0);

            const double ct = uMinus1_1(rng);
            const double st = std::sqrt(std::max(0.0, 1.0 - ct * ct));
            const double phi = uAngle(rng);
            const double x = r * st * std::cos(phi);
            const double y = r * st * std::sin(phi);
            const double z = r * ct;

            // Aarseth-Henon-Wielen 1974: rejection sampling по q = v / v_esc.
            const double v_esc = std::sqrt(2.0 * GM / std::sqrt(r * r + a2));
            double q;
            while (true) {
                q = u01(rng);
                const double g = q * q * std::pow(1.0 - q * q, 3.5);
                if (0.1 * u01(rng) < g) break;
            }
            const double v = q * v_esc;

            const double ct2 = uMinus1_1(rng);
            const double st2 = std::sqrt(std::max(0.0, 1.0 - ct2 * ct2));
            const double psi = uAngle(rng);
            const double vx = v * st2 * std::cos(psi);
            const double vy = v * st2 * std::sin(psi);
            const double vz = v * ct2;

            Object o(glm::dvec3(x, y, z), glm::dvec3(vx, vy, vz),
                     m_star, 5000.0, std::nullopt);
            o.Initalizing = false;
            o.radius = CalcRadiusAU(o.mass);
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
    mgr->registerScenario(std::make_unique<M13PlummerScenario>());
    return mgr;
}
