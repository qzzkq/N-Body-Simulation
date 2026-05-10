// Интеграционные тесты алгоритмов — ALG-01 по ALG-09
// Сравнивает BruteForce, Barnes-Hut CPU и Barnes-Hut CUDA.

#include <gtest/gtest.h>
#include <cmath>
#include <random>
#include <vector>
#include <glm/glm.hpp>

#include "object.hpp"
#include "bodysystem.hpp"
#include "physics.hpp"
#include "brut_force.hpp"
#include "barnes_hut.hpp"
#include "generators.hpp"

#ifdef USE_CUDA
#include "barnes_hut_cuda.cuh"
#endif

// ─── вспомогательные функции ──────────────────────────────────────────────────

static std::vector<Object> makeRandomBodies(int n, unsigned seed,
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

// Возвращает максимальное отклонение позиций между двумя наборами тел.
static double maxPosDiff(const std::vector<Object>& a, const std::vector<Object>& b) {
    double maxD = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
        maxD = std::max(maxD, glm::length(a[i].position - b[i].position));
    return maxD;
}

// Создаёт стабильную конфигурацию: центральная масса + N тел на круговых орбитах.
static std::vector<Object> makeRingConfig(int n) {
    std::vector<Object> objs;
    objs.emplace_back(glm::dvec3(0, 0, 0), glm::dvec3(0, 0, 0),
                      1.0, 1408.0, std::nullopt);  // центральная масса 1 M☉
    for (int i = 1; i < n; ++i) {
        const double r     = 1.0 + (i - 1) * 9.0 / (n - 2);
        const double v     = std::sqrt(physics::G / r);
        const double angle = i * 2.0 * physics::PI / (n - 1);
        objs.emplace_back(
            glm::dvec3(r * std::cos(angle), r * std::sin(angle), 0.0),
            glm::dvec3(-v * std::sin(angle),  v * std::cos(angle), 0.0),
            1e-4, 1000.0, std::nullopt);
    }
    return objs;
}

// Фикстура: сохраняет и восстанавливает глобальные параметры физики.
class AlgTest : public ::testing::Test {
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

// ─── ALG-01 ──────────────────────────────────────────────────────────────────
// 50 тел на стабильных орбитах, 100 шагов, dt = 0.001, theta = 0.3.
// Максимальное расхождение позиций BruteForce vs BH-CPU < 1e-3 AU.
// Используется кольцо, а не случайные тела: хаотичные сближения дают
// экспоненциальное расхождение, не связанное с точностью алгоритма.
TEST_F(AlgTest, BruteForceVsBH_CPU_50Bodies_100Steps) {
    physics::setSofteningAU(0.01);
    physics::setBarnesHutTheta(0.3);

    auto objs_bf = makeRingConfig(50);
    auto objs_bh = objs_bf;  // идентичная копия начального состояния

    const double dt = 0.001;
    for (int s = 0; s < 100; ++s) {
        simulationStepBrutForceCPU(objs_bf, dt, false);
        simulationStepBarnesHutCPU(objs_bh, dt, false);
    }

    const double diff = maxPosDiff(objs_bf, objs_bh);
    EXPECT_LT(diff, 1e-3) << "Макс. расхождение BF vs BH: " << diff << " AU";
}

// ─── ALG-02 ──────────────────────────────────────────────────────────────────
// BH-CPU vs BH-CUDA: одинаковые начальные условия, 100 шагов, dt = 0.001.
// Компилируется только при наличии CUDA (флаг USE_CUDA).
#ifdef USE_CUDA
TEST_F(AlgTest, BarnesHutCPU_vs_CUDA_50Bodies_100Steps) {
    physics::setSofteningAU(0.01);
    physics::setBarnesHutTheta(0.5);

    auto objs_cpu = makeRingConfig(50);
    auto objs_gpu = objs_cpu;  // идентичная копия

    const size_t n  = objs_cpu.size();
    const double dt = 0.001;

    InitBarnesHutCUDA(n);

    for (int s = 0; s < 100; ++s) {
        simulationStepBarnesHutCPU(objs_cpu,  dt, false);
        simulationStepBarnesHutCUDA(objs_gpu, dt, false, /*forceSync=*/true);
    }

    CleanupBarnesHutCUDA();

    const double diff = maxPosDiff(objs_cpu, objs_gpu);
    EXPECT_LT(diff, 1e-3) << "Макс. расхождение BH-CPU vs BH-CUDA: " << diff << " AU";
}
#endif  // USE_CUDA

// ─── ALG-03 ──────────────────────────────────────────────────────────────────
// Barnes-Hut с theta = 0 раскрывает все узлы дерева → результат близок к BruteForce.
// Допуск 1e-5 AU, а не машинная точность: разный порядок суммирования сил
// (попарный vs обход дерева) даёт накопленную ошибку FP ~1e-7 AU за один шаг.
TEST_F(AlgTest, BarnesHut_Theta0_EqualsBruteForce) {
    physics::setSofteningAU(0.01);
    physics::setBarnesHutTheta(0.0);

    auto objs_bf = makeRandomBodies(20, 9999);
    auto objs_bh = objs_bf;

    const double dt = 0.001;
    simulationStepBrutForceCPU(objs_bf, dt, false);
    simulationStepBarnesHutCPU(objs_bh, dt, false);

    const double diff = maxPosDiff(objs_bf, objs_bh);
    // theta=0 заставляет BH раскрывать все узлы, но порядок суммирования
    // в дереве отличается от последовательного попарного — отсюда FP-погрешность.
    EXPECT_LT(diff, 1e-5) << "theta=0 BH vs BF: " << diff << " AU";
}

// ─── ALG-04 ──────────────────────────────────────────────────────────────────
// Barnes-Hut: сохранение энергии при theta = 0.5, 1000 шагов, dt = 0.001.
// ΔE/E₀ < 0.5%.
TEST_F(AlgTest, BarnesHut_EnergyConservation) {
    physics::setSofteningAU(0.02);
    physics::setBarnesHutTheta(0.5);

    auto objs = makeRandomBodies(10, 42);

    const double E0 = physics::calculateTotalEnergy(objs);
    ASSERT_NE(E0, 0.0);

    const double dt = 0.001;
    for (int s = 0; s < 1000; ++s)
        simulationStepBarnesHutCPU(objs, dt, false);

    const double E1    = physics::calculateTotalEnergy(objs);
    const double drift = std::abs((E1 - E0) / E0);
    EXPECT_LT(drift, 5e-3) << "Дрейф энергии BH = " << drift * 100.0 << "%";
}

// ─── ALG-05 ──────────────────────────────────────────────────────────────────
// Barnes-Hut с одним телом: ускорение = (0,0,0), нет краша.
TEST_F(AlgTest, BarnesHut_SingleBody_NoForce) {
    std::vector<Object> objs;
    objs.emplace_back(glm::dvec3(1, 2, 3), glm::dvec3(0.5, 0, 0),
                      1.0, 1000.0, std::nullopt);

    ASSERT_NO_THROW(simulationStepBarnesHutCPU(objs, 0.001, false));

    // Единственное тело — никаких других масс → ускорение ноль
    EXPECT_NEAR(glm::length(objs[0].acceleration), 0.0, 1e-15);
}

// ─── ALG-06 ──────────────────────────────────────────────────────────────────
// N = 2: Barnes-Hut не делает приближений → результат совпадает с BruteForce.
// Ожидается расхождение < 1e-10 AU.
TEST_F(AlgTest, BarnesHut_N2_EqualsBruteForce) {
    physics::setSofteningAU(0.01);
    physics::setBarnesHutTheta(0.5);  // theta не важен при N=2

    auto objs_bf = makeRandomBodies(2, 55);
    auto objs_bh = objs_bf;

    const double dt = 0.001;
    simulationStepBrutForceCPU(objs_bf, dt, false);
    simulationStepBarnesHutCPU(objs_bh, dt, false);

    const double diff = maxPosDiff(objs_bf, objs_bh);
    EXPECT_LT(diff, 1e-10) << "N=2 BF vs BH: " << diff << " AU";
}

// ─── ALG-07 ──────────────────────────────────────────────────────────────────
// pause = true: состояние тел не изменяется ни для одного из алгоритмов.
TEST_F(AlgTest, Pause_BodiesDoNotMove) {
    auto objs_bf = makeRandomBodies(5, 77);
    auto objs_bh = objs_bf;

    // Сохраняем исходные позиции
    std::vector<glm::dvec3> before_bf(objs_bf.size()), before_bh(objs_bh.size());
    for (size_t i = 0; i < objs_bf.size(); ++i) {
        before_bf[i] = objs_bf[i].position;
        before_bh[i] = objs_bh[i].position;
    }

    simulationStepBrutForceCPU(objs_bf, 0.001, /*pause=*/true);
    simulationStepBarnesHutCPU(objs_bh, 0.001, /*pause=*/true);

    for (size_t i = 0; i < objs_bf.size(); ++i) {
        EXPECT_EQ(objs_bf[i].position, before_bf[i])
            << "BF сдвинул тело " << i << " при pause=true";
        EXPECT_EQ(objs_bh[i].position, before_bh[i])
            << "BH сдвинул тело " << i << " при pause=true";
    }
}

// ─── ALG-08 ──────────────────────────────────────────────────────────────────
// ScenarioManager: одинаковый seed порождает побитово идентичные системы.
TEST(ALG, ScenarioManager_SeedReproducibility) {
    auto mgr = CreateDefaultManager();
    ASSERT_TRUE(mgr->isValidIndex(0));

    GenParams params;
    params.count       = 50;
    params.centralMass = 1.0;
    params.baseMass    = 1e-3;
    params.minRadius   = 1.0f;
    params.maxRadius   = 5.0f;
    params.spread      = 0.1f;
    params.seed        = 42;

    std::vector<Object> run1, run2;
    mgr->runScenario(0, run1, params);
    mgr->runScenario(0, run2, params);

    ASSERT_EQ(run1.size(), run2.size());
    for (size_t i = 0; i < run1.size(); ++i) {
        EXPECT_EQ(run1[i].position, run2[i].position)
            << "Расхождение позиций у тела " << i;
        EXPECT_EQ(run1[i].velocity, run2[i].velocity)
            << "Расхождение скоростей у тела " << i;
        EXPECT_EQ(run1[i].mass, run2[i].mass)
            << "Расхождение масс у тела " << i;
    }
}

// ─── ALG-09 ──────────────────────────────────────────────────────────────────
// Кольцо из 100 тел (seed=42), симуляция 1 год.
// Все тела остаются в пределах 2×maxRadius от центра.
TEST_F(AlgTest, Ring_OrbitalStability_1Year) {
    physics::setSofteningAU(0.01);
    physics::setBarnesHutTheta(0.5);

    auto mgr = CreateDefaultManager();
    ASSERT_TRUE(mgr->isValidIndex(0));

    GenParams params;
    params.count       = 99;    // + 1 центральное тело = 100 итого
    params.centralMass = 1.0;
    params.baseMass    = 1e-4;
    params.minRadius   = 1.0f;
    params.maxRadius   = 5.0f;
    params.spread      = 0.1f;
    params.seed        = 42;

    std::vector<Object> objs;
    mgr->runScenario(0, objs, params);
    ASSERT_FALSE(objs.empty());

    BodySystem sys(objs);
    sys.transPointToSystem(objs);

    const double dt       = 1e-3;
    const int    steps    = static_cast<int>(1.0 / dt);
    const double maxBound = 2.0 * params.maxRadius;

    for (int s = 0; s < steps; ++s)
        simulationStepBarnesHutCPU(objs, dt, false);

    for (size_t i = 0; i < objs.size(); ++i) {
        const double r = glm::length(objs[i].position);
        EXPECT_LT(r, maxBound)
            << "Тело " << i << " покинуло систему: r = " << r << " AU";
    }
}
