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

    void generate(std::vector<Object>& out, const GenParams& p,
                  std::vector<GraphicState>* /*outGraphics*/) override {
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

    void generate(std::vector<Object>& out, const GenParams& p,
                  std::vector<GraphicState>* /*outGraphics*/) override {
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

    void generate(std::vector<Object>& out, const GenParams& p,
                  std::vector<GraphicState>* /*outGraphics*/) override {
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

class GalaxyCollisionScenario : public IScenario {
public:
    std::string getName() const override {
        return "Столкновение 2 галактик (синяя+красная)";
    }

    void generate(std::vector<Object>& out, const GenParams& p,
                  std::vector<GraphicState>* g) override {
        out.clear();
        if (g) g->clear();

        const int total = std::max(p.count, 0);
        const int n1 = total / 2;
        const int n2 = total - n1;
        out.reserve(static_cast<size_t>(total) + 2);
        if (g) g->reserve(static_cast<size_t>(total) + 2);

        std::mt19937 rng(p.seed ? p.seed : 1337u);

        const float R_max = std::max(p.maxRadius, 1.0f);
        const float D = R_max * 4.0f;
        const float impact = R_max * 0.5f;

        const double totalMass = 2.0 * p.centralMass + total * p.baseMass;
        const double vAppr = std::sqrt(physics::G * totalMass / (2.0 * D)) * 1.05;

        glm::vec3 c1(-D, 0.0f, +impact * 0.5f);
        glm::vec3 c2(+D, 0.0f, -impact * 0.5f);
        glm::vec3 vb1(+static_cast<float>(vAppr), 0.0f, 0.0f);
        glm::vec3 vb2(-static_cast<float>(vAppr), 0.0f, 0.0f);

        glm::mat3 R1(1.0f);
        const float ang2 = 0.6f;
        const float ca = std::cos(ang2), sa = std::sin(ang2);
        glm::mat3 R2(
            glm::vec3(1.0f, 0.0f, 0.0f),
            glm::vec3(0.0f,   ca,  -sa),
            glm::vec3(0.0f,   sa,   ca)
        );

        addGalaxy(out, g, rng, n1, p, c1, vb1, R1, /*hueBase=*/220.0f);
        addGalaxy(out, g, rng, n2, p, c2, vb2, R2, /*hueBase=*/  5.0f);
    }

private:
    static glm::vec3 hsv2rgb(float h, float s, float v) {
        h = std::fmod(h, 360.0f); if (h < 0) h += 360.0f;
        float c = v * s;
        float x = c * (1.0f - std::fabs(std::fmod(h / 60.0f, 2.0f) - 1.0f));
        float m = v - c;
        glm::vec3 rgb;
        if      (h <  60.0f) rgb = glm::vec3(c, x, 0);
        else if (h < 120.0f) rgb = glm::vec3(x, c, 0);
        else if (h < 180.0f) rgb = glm::vec3(0, c, x);
        else if (h < 240.0f) rgb = glm::vec3(0, x, c);
        else if (h < 300.0f) rgb = glm::vec3(x, 0, c);
        else                 rgb = glm::vec3(c, 0, x);
        return rgb + glm::vec3(m);
    }

    static void addGalaxy(std::vector<Object>& out,
                          std::vector<GraphicState>* g,
                          std::mt19937& rng,
                          int n,
                          const GenParams& p,
                          const glm::vec3& center,
                          const glm::vec3& bulkVel,
                          const glm::mat3& R,
                          float hueBase)
    {
        std::uniform_real_distribution<float> u01(0.0f, 1.0f);
        std::normal_distribution<float>       gauss(0.0f, 1.0f);
        std::uniform_real_distribution<float> uAngle(0.0f, 6.28318530718f);

        {
            Object bh(glm::dvec3(center), glm::dvec3(bulkVel), p.centralMass, 141000.0, std::nullopt);
            bh.Initalizing = false;
            bh.radius = CalcRadiusAU(bh.mass);
            out.push_back(bh);
            if (g) {
                GraphicState gs;
                glm::vec3 core = hsv2rgb(hueBase, 0.25f, 1.0f);
                gs.color = glm::vec4(glm::mix(core, glm::vec3(1.0f), 0.6f), 1.0f);
                g->push_back(gs);
            }
        }

        const float R_min = std::max(p.minRadius, 0.1f);
        const float R_max = std::max(p.maxRadius, R_min + 1.0f);
        const float R_d   = (R_min + R_max) * 0.25f;
        const float bulgeFrac   = 0.15f;
        const float bulgeRadius = std::max(R_min * 2.0f, R_max * 0.08f);

        for (int i = 0; i < n; ++i) {
            bool isBulge = u01(rng) < bulgeFrac;
            glm::vec3 pos(0.0f), vel(0.0f);

            if (isBulge) {
                // Изотропная сфера, плотность ~ 1/(1+r/a)^2
                float u = u01(rng);
                float r = bulgeRadius * u / std::max(1e-4f, 1.0f - u);
                r = std::min(r, bulgeRadius * 4.0f);
                float cosT = 1.0f - 2.0f * u01(rng);
                float sinT = std::sqrt(std::max(0.0f, 1.0f - cosT * cosT));
                float phi  = uAngle(rng);
                pos = glm::vec3(r * sinT * std::cos(phi),
                                r * cosT,
                                r * sinT * std::sin(phi));
                double dist = std::max(static_cast<double>(r), 1e-6);
                double vc = std::sqrt(physics::G * p.centralMass / dist);
                glm::vec3 randv(gauss(rng), gauss(rng), gauss(rng));
                glm::vec3 cross = glm::cross(pos, randv);
                float len = glm::length(cross);
                glm::vec3 tdir = (len > 1e-6f) ? (cross / len) : glm::vec3(0, 1, 0);
                vel = tdir * static_cast<float>(vc) * 0.6f;
            } else {
                // Экспоненциальный диск (Rayleigh-аппроксимация)
                float r;
                int guard = 0;
                do {
                    float u = std::max(u01(rng), 1e-6f);
                    r = R_d * std::sqrt(-2.0f * std::log(u));
                    if (++guard > 32) { r = std::clamp(r, R_min, R_max); break; }
                } while (r < R_min || r > R_max);

                float a = uAngle(rng);
                float zSigma = p.spread * std::exp(-r / (R_d * 2.5f));
                float z = gauss(rng) * std::max(zSigma, 0.05f);

                pos = glm::vec3(r * std::cos(a), z, r * std::sin(a));
                glm::vec3 tdir(-std::sin(a), 0.0f, std::cos(a));
                double dist = std::sqrt(static_cast<double>(r) * r + static_cast<double>(z) * z);
                double vc = std::sqrt(physics::G * p.centralMass / std::max(dist, 1e-9));
                vel = tdir * static_cast<float>(vc);
            }

            pos = R * pos;
            vel = R * vel;
            glm::dvec3 wp = glm::dvec3(pos) + glm::dvec3(center);
            glm::dvec3 wv = glm::dvec3(vel) + glm::dvec3(bulkVel);

            Object o(wp, wv, p.baseMass, 1410.0, std::nullopt);
            o.Initalizing = false;
            o.radius = CalcRadiusAU(o.mass);
            out.push_back(std::move(o));

            if (g) {
                GraphicState gs;
                float hue = hueBase + (u01(rng) * 2.0f - 1.0f) * 22.0f;
                float sat = 0.55f + u01(rng) * 0.40f;
                float val = 0.65f + u01(rng) * 0.35f;
                if (isBulge) {
                    val = std::min(1.0f, val + 0.15f);
                    sat = std::max(0.20f, sat - 0.25f);
                }
                glm::vec3 col = hsv2rgb(hue, sat, val);
                gs.color = glm::vec4(col, 1.0f);
                g->push_back(gs);
            }
        }
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

void ScenarioManager::runScenario(size_t index, std::vector<Object>& out, const GenParams& params,
                                  std::vector<GraphicState>* outGraphics) {
    if (isValidIndex(index)) {
        std::cout << "Генерируем: " << scenarios_[index]->getName() << "...\n";
        scenarios_[index]->generate(out, params, outGraphics);
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
    mgr->registerScenario(std::make_unique<GalaxyCollisionScenario>());
    return mgr;
}
