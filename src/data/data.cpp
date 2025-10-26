#include <vector>
#include <string>
#include <iostream>
#include <glm/glm.hpp>
#include "H5Cpp.h"

using namespace H5;

struct Particle {
    glm::dvec3 position;
    glm::dvec3 velocity;
    double     mass;
    double     radius;
};

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

