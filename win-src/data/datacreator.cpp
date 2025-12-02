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

    const double MASS_REF = 1.0e22;
    const double DIST_REF = 60.0;
    const double V_REF    = 1.0;

    const double MASS_HEAVY = MASS_REF;
    const double MASS_OUTER = MASS_REF * 0.001;

    const double R_HEAVY = 2.5;
    const double R_OUTER = 1.5;

    auto estimate_circular_velocity = [&](double r) {
        return V_REF * std::sqrt(2.0 * DIST_REF / r);
    };

    const double r_pair  = 20.0;
    const double R_orbit = r_pair * 0.5;

    double v_single = estimate_circular_velocity(r_pair);
    double v_binary = 0.5 * v_single;

    {
        Particle a;
        a.position = glm::dvec3(-R_orbit, 0.0, 0.0);
        a.velocity = glm::dvec3(0.0, -v_binary, 0.0);
        a.mass     = MASS_HEAVY;
        a.radius   = R_HEAVY;
        parts.push_back(a);
    }

    {
        Particle b;
        b.position = glm::dvec3(+R_orbit, 0.0, 0.0);
        b.velocity = glm::dvec3(0.0, +v_binary, 0.0);
        b.mass     = MASS_HEAVY;
        b.radius   = R_HEAVY;
        parts.push_back(b);
    }

    const double R_outer = 300.0;

    Particle c;
    c.position = glm::dvec3(0.0, R_outer, 0.0);

    double v_circ_base  = estimate_circular_velocity(R_outer);
    double v_circ_outer = std::sqrt(2.0) * v_circ_base;

    double v_outer = 1.1 * v_circ_outer;

    c.velocity = glm::dvec3(-v_outer, 0.0, 0.0);

    c.mass   = MASS_OUTER;
    c.radius = R_OUTER;
    parts.push_back(c);

    Writer("data/three_body_hierarchical_testparticle.h5", "Particles", parts);
    std::cout << "Saved " << parts.size()
              << " bodies to data/three_body_hierarchical_testparticle.h5\n";

    return 0;
}


