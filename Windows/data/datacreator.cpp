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

    const double initMass = 5.0 * std::pow(10.0, 20.0) / 5.0;
    const double rho_center = 141000.0;
    const double rho_sat    = 14100.0;

    std::vector<Particle> parts;
    parts.reserve(4);

    // центр
    // спутник 1
    {
        double m = initMass * 100.0;
        parts.push_back({
            { 60.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            m,
            radiusFromMassDensity(m, rho_sat)
        });
    }
    // спутник 2
    {
        double m = initMass * 100.0;
        parts.push_back({
            {0.0, 0.0, 0.0},
            {  0.0,-1.0, 0.0},
            m,
            radiusFromMassDensity(m, rho_sat)
        });
    }


    Writer("data/3bodies.h5", "Particles", parts);
    std::cout << "Saved " << parts.size() << " bodies to data/3bodies.h5\n";
    return 0;
}
