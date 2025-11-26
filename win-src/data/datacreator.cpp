#include <vector>
#include <string>
#include <iostream>
#include <glm/glm.hpp>
#include <filesystem>
#include <cmath>
#include <algorithm>
#include "H5Cpp.h"
#include "object.hpp"
using namespace H5;

struct Particle {
    glm::dvec3 position;
    glm::dvec3 velocity;
    double     mass;
    double     radius;
};

std::vector<Object> objs = {};
double initMass = 5.0f * std::pow(10.0f, 20.0f) / 5.0f;

static inline double radiusFromMassDensity(double m, double rho) {
    constexpr double pi = 3.14159265358979323846;
    const double r_m = std::cbrt((3.0 * m) / (4.0 * pi * rho));
    return r_m / 100000.0;
}


static CompType MakeFileType() {
    hsize_t v3[1] = {3};
    ArrayType vec3T(PredType::NATIVE_DOUBLE, 1, v3);

    constexpr size_t D    = sizeof(double);
    constexpr size_t POS  = 0;          
    constexpr size_t VEL  = POS + 3*D;  
    constexpr size_t MASS = VEL + 3*D;  
    constexpr size_t RAD  = MASS + D;   
    constexpr size_t SIZE = RAD + D;    

    CompType t(SIZE);
    t.insertMember("position", POS,  vec3T);
    t.insertMember("velocity", VEL,  vec3T);
    t.insertMember("mass",     MASS, PredType::NATIVE_DOUBLE);
    t.insertMember("radius",   RAD,  PredType::NATIVE_DOUBLE);
    return t;
}

static CompType MakeMemType() {
    hsize_t v3[1] = {3};
    ArrayType vec3T(PredType::NATIVE_DOUBLE, 1, v3);

    CompType t(sizeof(Particle));
    t.insertMember("position", HOFFSET(Particle, position), vec3T);
    t.insertMember("velocity", HOFFSET(Particle, velocity), vec3T);
    t.insertMember("mass",     HOFFSET(Particle, mass),     PredType::NATIVE_DOUBLE);
    t.insertMember("radius",   HOFFSET(Particle, radius),   PredType::NATIVE_DOUBLE);
    return t;
}

void Writer(const std::string& fileName,
            const std::string& dsetName,
            const std::vector<Particle>& particles)
{
    H5File file(fileName, H5F_ACC_TRUNC);
    CompType fileType = MakeFileType();
    CompType memType  = MakeMemType();

    hsize_t dims[1] = { static_cast<hsize_t>(particles.size()) };
    DataSpace space(1, dims);

    DataSet dset = file.createDataSet(dsetName, fileType, space);
    dset.write(particles.data(), memType);
}

std::vector<Particle> Reader(const std::string& fileName,
                             const std::string& dsetName)
{
    H5File file(fileName, H5F_ACC_RDONLY);
    DataSet dset = file.openDataSet(dsetName);

    DataSpace fspace = dset.getSpace();
    hsize_t dims[H5S_MAX_RANK]{};
    fspace.getSimpleExtentDims(dims);
    const hsize_t n = dims[0];

    std::vector<Particle> out(static_cast<size_t>(n));

    hsize_t mdims[1] = { n };
    DataSpace mspace(1, mdims);

    CompType memType = MakeMemType();
    dset.read(out.data(), memType, mspace, fspace);
    return out;
}

std::vector<std::string> ListH5Files(const std::string& dir) {
    std::vector<std::string> out;
    namespace fs = std::filesystem;
    if (!fs::exists(dir)) return out;
    for (const auto& e : fs::directory_iterator(dir)) {
        if (e.is_regular_file()) {
            const auto& p = e.path();
            if (p.has_extension() && p.extension() == ".h5")
                out.push_back(p.string());
        }
    }
    std::sort(out.begin(), out.end());
    return out;
}


int main() {
    std::vector<Particle> parts;

    // --- Реальные массы (кг), только для ОТНОСИТЕЛЬНЫХ соотношений ---
    const double mass_sun     = 1.9885e30;
    const double mass_mercury = 3.3011e23;
    const double mass_venus   = 4.8675e24;
    const double mass_earth   = 5.97237e24;
    const double mass_mars    = 6.4171e23;
    const double mass_jupiter = 1.8982e27;
    const double mass_saturn  = 5.6834e26;
    const double mass_uranus  = 8.6810e25;
    const double mass_neptune = 1.02413e26;

    // Средние расстояния от Солнца (в а.е.)
    const double a_mercury = 0.387;
    const double a_venus   = 0.723;
    const double a_earth   = 1.0;
    const double a_mars    = 1.524;
    const double a_jupiter = 5.203;
    const double a_saturn  = 9.537;
    const double a_uranus  = 19.191;
    const double a_neptune = 30.07;

    // --- Эталонные единицы из рабочего 3bodies-примера ---
    const double MASS_REF = 1.0e22;  // масса "типичного" тела в исходном примере
    const double DIST_REF = 60.0;    // расстояние между телами в исходном примере
    const double V_REF    = 1.0;     // скорость в исходном примере

    // Масштаб масс: делаем так, чтобы Солнце имело массу MASS_REF
    const double MASS_SCALE = MASS_REF / mass_sun;

    auto mass_sim = [&](double m_real) {
        return m_real * MASS_SCALE;
    };

    // Масштаб расстояния: Земля -> радиус = DIST_REF (как в 3bodies)
    const double DIST_SCALE = DIST_REF / a_earth;  // = 60.0

    // Радиусы чисто визуальные
    const double radius_sun_sim    = 4.0;
    const double radius_planet_sim = 1.0;

    // Коэффициент по скорости: чуть меньше "круговой" — устойчивые эллипсы
    const double V_FACTOR = 0.8;

    // Приближённая оценка "круговой" скорости в тех же юнитах, что и в 3bodies:
    // v_circ ≈ V_REF * sqrt(2 * DIST_REF / r)
    auto estimate_circular_velocity = [&](double r_sim) {
        return V_REF * std::sqrt(2.0 * DIST_REF / r_sim);
    };

    // --- Солнце в центре ---
    {
        Particle sun;
        sun.position = glm::dvec3(0.0, 0.0, 0.0);
        sun.velocity = glm::dvec3(0.0, 0.0, 0.0);
        sun.mass     = mass_sim(mass_sun);  // == MASS_REF
        sun.radius   = radius_sun_sim;
        parts.push_back(sun);
    }

    // Хелпер для добавления планет на (почти) круговые орбиты
    auto add_planet = [&](double a_AU, double m_real) {
        double r_sim  = a_AU * DIST_SCALE;          // расстояние в сим-юнитах
        double v_circ = estimate_circular_velocity(r_sim);
        double v_sim  = V_FACTOR * v_circ;          // чуть меньше круговой

        Particle p;
        p.position = glm::dvec3(r_sim, 0.0, 0.0);   // на оси X
        p.velocity = glm::dvec3(0.0,  v_sim, 0.0);  // скорость вдоль +Y
        p.mass     = mass_sim(m_real);
        p.radius   = radius_planet_sim;
        parts.push_back(p);
    };

    // --- Планеты ---
    add_planet(a_mercury, mass_mercury);
    add_planet(a_venus,   mass_venus);
    add_planet(a_earth,   mass_earth);
    add_planet(a_mars,    mass_mars);
    add_planet(a_jupiter, mass_jupiter);
    add_planet(a_saturn,  mass_saturn);
    add_planet(a_uranus,  mass_uranus);
    add_planet(a_neptune, mass_neptune);

    // --- ТЯЖЁЛОЕ ПРОЛЁТНОЕ ТЕЛО ("чёрная дыра" или звезда-нарушитель) ---

    // Пусть оно в 5 раз тяжелее Солнца по массе:
    const double mass_intruder_real = 5.0 * mass_sun;

    Particle intr;
    intr.mass   = mass_sim(mass_intruder_real); // ~ 5 * MASS_REF
    intr.radius = 3.0;

    // Стартовая позиция: снаружи орбиты Сатурна, на окраине системы
    // Немного смещаем по X и Y, чтобы траектория прошла "насквозь"
    intr.position = glm::dvec3(-500.0, 300.0, 0.0);

    // Скорость направлена примерно к центру,
    // по величине сравнима со скоростями планет.
    intr.velocity = glm::dvec3(1.0, -0.5, 0.0);

    parts.push_back(intr);

    // --- Запись HDF5 ---
    Writer("data/solar_system_intruder.h5", "Particles", parts);
    std::cout << "Saved " << parts.size()
              << " bodies to data/solar_system_intruder.h5\n";

    return 0;
}






