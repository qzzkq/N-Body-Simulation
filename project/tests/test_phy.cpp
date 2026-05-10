// Юнит-тесты физики — PHY-01 по PHY-08
// Единицы измерения: AU, M☉, юлианский год.  G = 4π² ≈ 39.478 AU³/(M☉·год²).

#include <gtest/gtest.h>
#include <cmath>
#include <random>
#include <vector>
#include <glm/glm.hpp>

#include "object.hpp"
#include "bodysystem.hpp"
#include "physics.hpp"
#include "brut_force.hpp"

// ─── вспомогательные функции ──────────────────────────────────────────────────

// Создаёт N тел со случайными позициями, скоростями и массами с фиксированным seed.
static std::vector<Object> makeRandomBodies(int n, unsigned seed = 42,
                                             double posRange = 5.0,
                                             double velRange = 1.0) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> pd(-posRange, posRange);
    std::uniform_real_distribution<double> vd(-velRange, velRange);
    std::uniform_real_distribution<double> md(0.5, 2.0);

    std::vector<Object> objs;
    objs.reserve(n);
    for (int i = 0; i < n; ++i)
        objs.emplace_back(
            glm::dvec3(pd(rng), pd(rng), pd(rng)),
            glm::dvec3(vd(rng), vd(rng), vd(rng)),
            md(rng), 1000.0, std::nullopt);
    return objs;
}

// Фикстура: сохраняет и восстанавливает глобальные параметры физики до/после каждого теста.
class PhysicsTest : public ::testing::Test {
protected:
    void SetUp() override {
        savedSoftening = physics::getSofteningAU();
        savedTheta     = physics::getBarnesHutTheta();
    }
    void TearDown() override {
        physics::setSofteningAU(savedSoftening);
        physics::setBarnesHutTheta(savedTheta);
    }
    double savedSoftening{}, savedTheta{};
};

// ─── PHY-01 ──────────────────────────────────────────────────────────────────
// Два тела, m = 1 M☉, расстояние = 1 AU.  Один шаг BruteForce, dt = 1e-6.
// Ожидается: |a| = G·m/r² = 4π² ≈ 39.478 AU/год²  (погрешность < 0.001%).
// Интегратор Yoshida4: за один шаг из состояния покоя v ≈ a·dt, поэтому a = v/dt.
TEST_F(PhysicsTest, GravAccel_TwoBodies_1AU) {
    physics::setSofteningAU(0.0);

    std::vector<Object> objs;
    objs.emplace_back(glm::dvec3(0, 0, 0), glm::dvec3(0, 0, 0),
                      1.0, 1408.0, std::nullopt);
    objs.emplace_back(glm::dvec3(1, 0, 0), glm::dvec3(0, 0, 0),
                      1.0, 1408.0, std::nullopt);

    const double dt = 1e-6;
    simulationStepBrutForceCPU(objs, dt, false);

    // Yoshida4 сохраняет |v| ≈ a·dt для первого шага (сумма коэффициентов = 1).
    const double accel = glm::length(objs[0].velocity) / dt;
    EXPECT_NEAR(accel, physics::G, physics::G * 1e-5);   // < 0.001%
}

// ─── PHY-02 ──────────────────────────────────────────────────────────────────
// Структурированная конфигурация: Солнце + 9 планет на круговых орбитах 1–9 AU.
// 1000 шагов, dt = 0.001.  ΔE/E₀ < 0.1%.
// Случайные тела не используются — хаотичные сближения искажают меру сохранения энергии.
TEST_F(PhysicsTest, EnergyConservation_BruteForce) {
    physics::setSofteningAU(0.0);

    // Центральное Солнце + 9 лёгких планет на точно круговых орбитах
    std::vector<Object> objs;
    objs.emplace_back(glm::dvec3(0, 0, 0), glm::dvec3(0, 0, 0),
                      1.0, 1408.0, std::nullopt);

    for (int i = 1; i <= 9; ++i) {
        const double r = i * 1.0;
        const double v = std::sqrt(physics::G / r);  // скорость круговой орбиты
        objs.emplace_back(
            glm::dvec3(r, 0, 0),
            glm::dvec3(0, v, 0),
            1e-4, 1000.0, std::nullopt);
    }

    BodySystem sys(objs);
    sys.transPointToSystem(objs);

    const double E0 = physics::calculateTotalEnergy(objs);
    ASSERT_NE(E0, 0.0) << "Начальная энергия не может быть равна нулю";

    const double dt = 0.001;
    for (int s = 0; s < 1000; ++s)
        simulationStepBrutForceCPU(objs, dt, false);

    const double E1    = physics::calculateTotalEnergy(objs);
    const double drift = std::abs((E1 - E0) / E0);
    EXPECT_LT(drift, 1e-3) << "Дрейф энергии = " << drift * 100.0 << "%";
}

// ─── PHY-03 ──────────────────────────────────────────────────────────────────
// После transPointToSystem суммарный импульс Σ(mᵢ·vᵢ) < 1e-12 M☉·AU/год.
TEST(PHY, MomentumConservation_AfterBoost) {
    auto objs = makeRandomBodies(8, 7);

    BodySystem sys(objs);
    sys.transPointToSystem(objs);

    glm::dvec3 p(0.0);
    for (const auto& o : objs) p += o.mass * o.velocity;

    EXPECT_NEAR(glm::length(p), 0.0, 1e-12);
}

// ─── PHY-04 ──────────────────────────────────────────────────────────────────
// Солнце: m = 1 M☉, ρ = 1408 кг/м³ → r ≈ 0.00465 AU  (погрешность < 1%).
TEST(PHY, CalculateRadius_Sun) {
    const float r = physics::calculateRadius(1.0, 1408.0);
    EXPECT_NEAR(r, 0.00465f, 0.00465f * 0.01f);
}

// ─── PHY-05 ──────────────────────────────────────────────────────────────────
// Нулевая или отрицательная плотность → возвращает 1.0f без краша.
TEST(PHY, CalculateRadius_ZeroOrNegDensity) {
    EXPECT_FLOAT_EQ(physics::calculateRadius(1.0,  0.0), 1.0f);
    EXPECT_FLOAT_EQ(physics::calculateRadius(1.0, -5.0), 1.0f);
}

// ─── PHY-06 ──────────────────────────────────────────────────────────────────
// Земля на r = 1 AU, v = 2π AU/год (круговая орбита вокруг 1 M☉).
// После ровно 1 года: отклонение позиции < 0.01 AU.
TEST_F(PhysicsTest, EarthOrbitalPeriod) {
    physics::setSofteningAU(0.0);

    std::vector<Object> objs;
    objs.emplace_back(glm::dvec3(0, 0, 0), glm::dvec3(0, 0, 0),
                      1.0, 1408.0, std::nullopt);                    // Солнце
    objs.emplace_back(glm::dvec3(1, 0, 0),
                      glm::dvec3(0, 2.0 * physics::PI, 0),
                      3e-6, 5000.0, std::nullopt);                   // Земля

    // Переход в систему центра масс (почти нет эффекта: масса Земли ≈ 0)
    BodySystem sys(objs);
    sys.transPointToSystem(objs);

    const glm::dvec3 pos0 = objs[1].position;

    const double dt    = 1e-4;
    const int    steps = static_cast<int>(1.0 / dt);  // 1 год
    for (int i = 0; i < steps; ++i)
        simulationStepBrutForceCPU(objs, dt, false);

    const double dev = glm::length(objs[1].position - pos0);
    EXPECT_LT(dev, 0.01) << "Отклонение орбиты за 1 год: " << dev << " AU";
}

// ─── PHY-07 ──────────────────────────────────────────────────────────────────
// Два тела при r = 0.001 AU: нет NaN/Inf при сглаживании eps=0 и eps=0.1.
TEST_F(PhysicsTest, Softening_NoSingularity) {
    for (double eps : {0.0, 0.1}) {
        physics::setSofteningAU(eps);

        std::vector<Object> objs;
        objs.emplace_back(glm::dvec3(0,     0, 0), glm::dvec3(0, 0, 0),
                          1.0, 1000.0, std::nullopt);
        objs.emplace_back(glm::dvec3(0.001, 0, 0), glm::dvec3(0, 0, 0),
                          1.0, 1000.0, std::nullopt);

        ASSERT_NO_THROW(simulationStepBrutForceCPU(objs, 1e-5, false))
            << "Исключение при eps=" << eps;

        for (const auto& o : objs) {
            EXPECT_FALSE(std::isnan(o.position.x)) << "NaN в позиции при eps=" << eps;
            EXPECT_FALSE(std::isinf(o.position.x)) << "Inf в позиции при eps=" << eps;
            EXPECT_FALSE(std::isnan(o.velocity.x)) << "NaN в скорости при eps=" << eps;
            EXPECT_FALSE(std::isinf(o.velocity.x)) << "Inf в скорости при eps=" << eps;
        }
    }
}

// ─── PHY-08 ──────────────────────────────────────────────────────────────────
// После transPointToSystem пересчёт BodySystem показывает getVel() ≈ (0,0,0).
TEST(PHY, BodySystem_VelZeroAfterBoost) {
    std::vector<Object> objs;
    objs.emplace_back(glm::dvec3( 2, 0, 0), glm::dvec3( 1,  0, 0), 2.0, 1000.0, std::nullopt);
    objs.emplace_back(glm::dvec3(-2, 0, 0), glm::dvec3(-1,  0, 0), 2.0, 1000.0, std::nullopt);
    objs.emplace_back(glm::dvec3( 0, 3, 0), glm::dvec3( 0,  1, 0), 1.0, 1000.0, std::nullopt);

    BodySystem sys(objs);
    sys.transPointToSystem(objs);

    // Пересчитываем параметры после буста — скорость ЦМ должна быть нулевой
    BodySystem sys2(objs);
    EXPECT_NEAR(glm::length(sys2.getVel()), 0.0, 1e-12);

    // Суммарный импульс тоже обнулился
    glm::dvec3 p(0.0);
    double M = 0.0;
    for (const auto& o : objs) { p += o.mass * o.velocity; M += o.mass; }
    EXPECT_NEAR(glm::length(p), 0.0, 1e-12);
}
