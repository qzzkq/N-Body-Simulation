// Юнит-тесты ввода/вывода — IO-01 по IO-10

#include <gtest/gtest.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <glm/glm.hpp>

#include "H5Cpp.h"
#include "object.hpp"
#include "graphic_state.hpp"
#include "data.hpp"
#include "physics.hpp"

namespace fs = std::filesystem;

// ─── фикстура: временная директория ──────────────────────────────────────────

class IOTest : public ::testing::Test {
protected:
    fs::path tmpDir;

    void SetUp() override {
        char buf[] = "/tmp/nbody_io_XXXXXX";
        ASSERT_NE(mkdtemp(buf), nullptr) << "mkdtemp провалился";
        tmpDir = buf;
    }

    void TearDown() override {
        fs::remove_all(tmpDir);
    }

    // Строит полный путь к файлу внутри временной директории.
    std::string p(const std::string& name) const {
        return (tmpDir / name).string();
    }

    // Создаёт набор тестовых тел с заданными позициями и массами.
    static std::vector<Object> makeObjects(int n = 5) {
        std::vector<Object> objs;
        for (int i = 0; i < n; ++i)
            objs.emplace_back(
                glm::dvec3(i * 1.5, i * 0.3, 0.0),
                glm::dvec3(0.0, 0.1 * i, 0.0),
                1.0 + 0.2 * i, 1000.0, std::nullopt);
        return objs;
    }

    // Создаёт набор Particle для тестов Writer/Reader.
    static std::vector<Particle> makeParticles(int n = 5) {
        std::vector<Particle> particles(n);
        for (int i = 0; i < n; ++i) {
            particles[i].position = glm::dvec3(i * 1.5, i * 0.3, 0.0);
            particles[i].velocity = glm::dvec3(0.0, 0.1 * i, 0.0);
            particles[i].mass     = 1.0 + 0.2 * i;
            particles[i].radius   = 0.01f;
            particles[i].color    = glm::vec4(0.8f, 0.2f, 0.1f + 0.05f * i, 1.0f);
        }
        return particles;
    }
};

// ─── IO-01 ───────────────────────────────────────────────────────────────────
// Writer + Reader: позиции сохраняются и восстанавливаются с точностью < 1e-5 AU.
TEST_F(IOTest, Writer_Reader_PositionRoundTrip) {
    const auto particles = makeParticles(10);
    const std::string fname = p("particles.h5");

    Writer(fname, "Particles", particles);

    const auto loaded = Reader(fname, "Particles");
    ASSERT_EQ(loaded.size(), particles.size());

    for (size_t i = 0; i < particles.size(); ++i) {
        EXPECT_NEAR(loaded[i].position.x, particles[i].position.x, 1e-5)
            << "Расхождение x у тела " << i;
        EXPECT_NEAR(loaded[i].position.y, particles[i].position.y, 1e-5)
            << "Расхождение y у тела " << i;
        EXPECT_NEAR(loaded[i].position.z, particles[i].position.z, 1e-5)
            << "Расхождение z у тела " << i;
    }
}

// ─── IO-02 ───────────────────────────────────────────────────────────────────
// CreateSimulationFile: атрибуты num_bodies и dt записаны корректно.
TEST_F(IOTest, CreateSimulationFile_Attributes) {
    const std::string fname  = p("sim.h5");
    const std::size_t nBodies = 7;
    const double      dt      = 0.00137;

    auto file = CreateSimulationFile(fname, nBodies, dt);
    CloseSimulationFile(file, 0);

    // Читаем атрибуты напрямую через H5Cpp для проверки
    H5::H5File f(fname, H5F_ACC_RDONLY);

    hsize_t nb = 0;
    f.openAttribute("num_bodies").read(H5::PredType::NATIVE_HSIZE, &nb);
    EXPECT_EQ(nb, nBodies);

    double dt_read = 0.0;
    f.openAttribute("dt").read(H5::PredType::NATIVE_DOUBLE, &dt_read);
    EXPECT_NEAR(dt_read, dt, 1e-15);
}

// ─── IO-03 ───────────────────────────────────────────────────────────────────
// CloseSimulationFile: атрибут num_frames равен переданному значению.
TEST_F(IOTest, CloseSimulationFile_NumFrames) {
    auto objs = makeObjects(4);
    std::vector<GraphicState> gfx(4);
    const std::string fname = p("frames.h5");

    auto file = CreateSimulationFile(fname, objs.size(), 0.001);
    for (std::size_t fr = 0; fr < 5; ++fr)
        WriteSimulationFrame(file, objs, gfx, fr);
    CloseSimulationFile(file, 5);

    H5::H5File f(fname, H5F_ACC_RDONLY);
    hsize_t nf = 0;
    f.openAttribute("num_frames").read(H5::PredType::NATIVE_HSIZE, &nf);
    EXPECT_EQ(nf, static_cast<hsize_t>(5));
}

// ─── IO-04 ───────────────────────────────────────────────────────────────────
// LoadObjectsFromFile: масса тел восстанавливается с точностью < 1e-10 M☉.
TEST_F(IOTest, LoadObjectsFromFile_MassRoundTrip) {
    const auto particles = makeParticles(3);
    const std::string fname = p("mass_rt.h5");

    Writer(fname, "Particles", particles);

    std::vector<Object> outObjs;
    ASSERT_TRUE(LoadObjectsFromFile(fname, "Particles", outObjs));
    ASSERT_EQ(outObjs.size(), particles.size());

    for (size_t i = 0; i < particles.size(); ++i)
        EXPECT_NEAR(outObjs[i].mass, particles[i].mass, 1e-10)
            << "Расхождение масс у тела " << i;
}

// ─── IO-05 ───────────────────────────────────────────────────────────────────
// LoadSystemFromTextFile: все поля разбираются корректно.
// Реальный формат файла: name massKg density_kg_m3 px_m py_m pz_m vx_ms vy_ms vz_ms cr cg cb
// Позиции в м → AU, скорости в м/с → AU/год, масса в кг → M☉ выполняются внутри загрузчика.
TEST_F(IOTest, LoadSystemFromTextFile_CorrectParsing) {
    using namespace physics;

    const std::string fname = p("system.txt");

    // Пересчитываем SI-единицы из астрономических, как это делает SaveSystemToTextFile
    const double sun_kg    = 1.0 / MASS_TO_SOLAR;
    const double earth_kg  = 3e-6 / MASS_TO_SOLAR;
    const double earth_m   = 1.0 / METERS_TO_AU;           // 1 AU в метрах
    const double earth_vms = 2.0 * PI / VELOCITY_TO_AU_PER_YEAR; // 2π AU/год → м/с
    const double mars_kg   = 3.2e-7 / MASS_TO_SOLAR;
    const double mars_m    = 1.52 / METERS_TO_AU;

    {
        std::ofstream f(fname);
        f << std::scientific << std::setprecision(15);
        f << "Sun   " << sun_kg   << " 1408.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 0.8\n";
        f << "Earth " << earth_kg << " 5000.0 " << earth_m << " 0.0 0.0 0.0 "
          << earth_vms << " 0.0 0.2 0.5 0.9\n";
        f << "Mars  " << mars_kg  << " 3900.0 " << mars_m << " 0.0 0.0 0.0 0.0 0.0 0.8 0.3 0.1\n";
    }

    std::vector<Object> objs;
    std::vector<GraphicState> gfx;
    ASSERT_TRUE(LoadSystemFromTextFile(fname, objs, &gfx));
    ASSERT_EQ(objs.size(), 3u);

    EXPECT_EQ(objs[0].name, "Sun");
    EXPECT_NEAR(objs[0].mass,    1.0,    1e-6);   // M☉
    EXPECT_NEAR(objs[0].density, 1408.0, 1e-3);

    EXPECT_EQ(objs[1].name, "Earth");
    EXPECT_NEAR(objs[1].position.x, 1.0,           1e-6);  // AU
    EXPECT_NEAR(objs[1].velocity.y, 2.0 * PI,      1e-6);  // AU/год
    EXPECT_NEAR(objs[1].mass,       3e-6,           1e-12); // M☉

    EXPECT_EQ(objs[2].name, "Mars");
    EXPECT_NEAR(objs[2].position.x, 1.52, 1e-6);
}

// ─── IO-06 ───────────────────────────────────────────────────────────────────
// Строки с комментариями (#) и пустые строки игнорируются.
TEST_F(IOTest, LoadSystemFromTextFile_CommentsIgnored) {
    using namespace physics;
    const std::string fname = p("comments.txt");
    const double sun_kg   = 1.0 / MASS_TO_SOLAR;
    const double earth_kg = 3e-6 / MASS_TO_SOLAR;
    const double earth_m  = 1.0 / METERS_TO_AU;
    const double earth_v  = 2.0 * PI / VELOCITY_TO_AU_PER_YEAR;
    {
        std::ofstream f(fname);
        f << std::scientific << std::setprecision(15);
        f << "# Это комментарий\n";
        f << "\n";
        f << "  \n";  // строка только из пробелов
        f << "Sun " << sun_kg << " 1408.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 0.8\n";
        f << "# Ещё один комментарий\n";
        f << "Earth " << earth_kg << " 5000.0 " << earth_m
          << " 0.0 0.0 0.0 " << earth_v << " 0.0 0.2 0.5 0.9\n";
    }

    std::vector<Object> objs;
    ASSERT_TRUE(LoadSystemFromTextFile(fname, objs));
    EXPECT_EQ(objs.size(), 2u);
    EXPECT_EQ(objs[0].name, "Sun");
    EXPECT_EQ(objs[1].name, "Earth");
}

// ─── IO-07 ───────────────────────────────────────────────────────────────────
// SaveSystemToTextFile + LoadSystemFromTextFile: обратное преобразование.
// Погрешность позиций < 1e-6 AU (потери на двойном преобразовании AU↔м).
TEST_F(IOTest, SaveLoadTextFile_RoundTrip) {
    auto origObjs = makeObjects(4);
    for (int i = 0; i < 4; ++i)
        origObjs[i].name = "Body" + std::to_string(i);

    const std::string fname = p("rt.txt");
    ASSERT_TRUE(SaveSystemToTextFile(fname, origObjs));

    std::vector<Object> loaded;
    ASSERT_TRUE(LoadSystemFromTextFile(fname, loaded));
    ASSERT_EQ(loaded.size(), origObjs.size());

    for (size_t i = 0; i < origObjs.size(); ++i) {
        EXPECT_NEAR(loaded[i].position.x, origObjs[i].position.x, 1e-6)
            << "Расхождение pos.x у тела " << i;
        EXPECT_NEAR(loaded[i].mass, origObjs[i].mass, 1e-10)
            << "Расхождение массы у тела " << i;
    }
}

// ─── IO-08 ───────────────────────────────────────────────────────────────────
// ListH5Files на пустой директории: возвращает пустой вектор без исключений.
TEST_F(IOTest, ListH5Files_EmptyDir) {
    const std::string emptyDir = p("empty");
    fs::create_directory(emptyDir);

    std::vector<std::string> files;
    EXPECT_NO_THROW(files = ListH5Files(emptyDir));
    EXPECT_TRUE(files.empty());
}

// ─── IO-09 ───────────────────────────────────────────────────────────────────
// LoadObjectsFromFile с несуществующим путём: возвращает false, outObjs пустой,
// исключение не пробрасывается (H5::Exception ловится внутри функции).
TEST_F(IOTest, LoadObjectsFromFile_NonexistentFile) {
    std::vector<Object> outObjs;
    bool result = true;
    EXPECT_NO_THROW(result = LoadObjectsFromFile(p("does_not_exist.h5"),
                                                  "Particles", outObjs));
    EXPECT_FALSE(result);
    EXPECT_TRUE(outObjs.empty());
}

// ─── IO-10 ───────────────────────────────────────────────────────────────────
// Цвета тел: сохраняем RGBA, загружаем обратно — значения совпадают.
TEST_F(IOTest, Colors_RoundTrip) {
    auto particles = makeParticles(4);
    const std::string fname = p("colors.h5");

    Writer(fname, "Particles", particles);

    bool hadColor = false;
    const auto loaded = Reader(fname, "Particles", &hadColor);

    EXPECT_TRUE(hadColor) << "Флаг наличия цветов должен быть true";
    ASSERT_EQ(loaded.size(), particles.size());

    for (size_t i = 0; i < particles.size(); ++i) {
        EXPECT_NEAR(loaded[i].color.r, particles[i].color.r, 1e-5f)
            << "Расхождение R у тела " << i;
        EXPECT_NEAR(loaded[i].color.g, particles[i].color.g, 1e-5f)
            << "Расхождение G у тела " << i;
        EXPECT_NEAR(loaded[i].color.b, particles[i].color.b, 1e-5f)
            << "Расхождение B у тела " << i;
    }
}
