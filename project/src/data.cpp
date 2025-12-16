#include "data.hpp"
#include "object.hpp"

#include <vector>
#include <string>
#include <iostream>
#include <filesystem>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <cstdio>

#include "H5Cpp.h"

using namespace H5;
namespace fs = std::filesystem;

static CompType MakeParticleType()
{
    hsize_t v3dims[1] = {3};
    ArrayType vec3Type(PredType::NATIVE_DOUBLE, 1, v3dims);

    CompType t(sizeof(Particle));
    t.insertMember("position", HOFFSET(Particle, position), vec3Type);
    t.insertMember("velocity", HOFFSET(Particle, velocity), vec3Type);
    t.insertMember("mass",     HOFFSET(Particle, mass),     PredType::NATIVE_DOUBLE);
    t.insertMember("radius",   HOFFSET(Particle, radius),   PredType::NATIVE_DOUBLE);
    return t;
}

static CompType& GetParticleType()
{
    static CompType t = MakeParticleType();
    return t;
}

void Writer(const std::string& fileName,
            const std::string& dsetName,
            const std::vector<Particle>& particles)
{
    try {
        H5File file(fileName, H5F_ACC_TRUNC);

        hsize_t dims[1] = { static_cast<hsize_t>(particles.size()) };
        DataSpace space(1, dims);

        CompType& type = GetParticleType();

        DataSet dset = file.createDataSet(dsetName, type, space);
        if (!particles.empty()) {
            dset.write(particles.data(), type);
        }

        file.flush(H5F_SCOPE_GLOBAL);
    }
    catch (const H5::Exception& e) {
        std::cerr << "HDF5 Writer error: " << e.getDetailMsg() << "\n";
        throw;
    }
}

std::vector<Particle> Reader(const std::string& fileName,
                             const std::string& dsetName)
{
    try {
        H5File file(fileName, H5F_ACC_RDONLY);
        DataSet dset = file.openDataSet(dsetName);

        DataSpace space = dset.getSpace();
        int rank = space.getSimpleExtentNdims();
        if (rank != 1) {
            throw std::runtime_error("Reader: dataset rank != 1");
        }

        hsize_t dims[1] = {0};
        space.getSimpleExtentDims(dims);
        std::size_t count = static_cast<std::size_t>(dims[0]);

        std::vector<Particle> result(count);
        if (count > 0) {
            CompType& type = GetParticleType();
            dset.read(result.data(), type);
        }

        return result;
    }
    catch (const H5::Exception& e) {
        std::cerr << "HDF5 Reader error: " << e.getDetailMsg() << "\n";
        throw;
    }
    catch (const std::exception& e) {
        std::cerr << "Reader error: " << e.what() << "\n";
        throw;
    }
}

std::vector<std::string> ListH5Files(const std::string& dir)
{
    std::vector<std::string> files;

    try {
        if (!fs::exists(dir)) {
            return files;
        }

        for (const auto& entry : fs::directory_iterator(dir)) {
            if (!entry.is_regular_file()) continue;
            const auto& path = entry.path();
            if (path.extension() == ".h5") {
                files.push_back(path.string());
            }
        }

        std::sort(files.begin(), files.end());
    }
    catch (const std::exception& e) {
        std::cerr << "ListH5Files error: " << e.what() << "\n";
    }

    return files;
}

bool LoadObjectsFromFile(const std::string& filePath,
                         const std::string& dsetName,
                         std::vector<Object>& outObjs)
{
    std::vector<Particle> parts;
    try {
        parts = Reader(filePath, dsetName);
    }
    catch (const std::exception& e) {
        std::cerr << "LoadObjectsFromFile: failed to read HDF5: "
                  << e.what() << "\n";
        return false;
    }

    if (parts.empty()) {
        std::cerr << "LoadObjectsFromFile: dataset is empty\n";
        return false;
    }

    outObjs.clear();
    outObjs.reserve(parts.size());

    for (const auto& p : parts) {
        Object o(p.position, p.velocity,
                 p.mass,            // mass
                 1410.0f,           // density (как в spawnSystem)
                 std::nullopt);     // radius – позже уточним

        o.Initalizing = false;

        if (std::isfinite(p.radius) && p.radius > 0.0) {
            o.radius = static_cast<float>(p.radius);
        } else {
            if (o.density > 0.0) {
                constexpr double pi = 3.14159265358979323846;
                const double r_m = std::cbrt((3.0 * o.mass) / (4.0 * pi * o.density));
                o.radius = static_cast<float>(r_m / 100000.0);
            } else {
                o.radius = 1.0f;
            }
        }
        outObjs.push_back(std::move(o));
    }

    return true;
}


struct InitRecord {
    glm::dvec3 position;
    double     mass;
    double     radius;
};

static CompType GetInitType() {
    hsize_t v3dims[1] = {3};
    ArrayType vec3Type(PredType::NATIVE_DOUBLE, 1, v3dims);
    CompType t(sizeof(InitRecord));
    t.insertMember("position", HOFFSET(InitRecord, position), vec3Type);
    t.insertMember("mass",     HOFFSET(InitRecord, mass),     PredType::NATIVE_DOUBLE);
    t.insertMember("radius",   HOFFSET(InitRecord, radius),   PredType::NATIVE_DOUBLE);
    return t;
}

static ArrayType GetPosType() {
    hsize_t v3dims[1] = {3};
    return ArrayType(PredType::NATIVE_DOUBLE, 1, v3dims);
}

H5::H5File CreateSimulationFile(const std::string& fileName, std::size_t numBodies, double dt)
{
    H5File file(fileName, H5F_ACC_TRUNC);

    hsize_t scalarDims[1] = {1};
    DataSpace scalarSpace(1, scalarDims);

    hsize_t nb = static_cast<hsize_t>(numBodies);
    Attribute attNb = file.createAttribute("num_bodies", PredType::NATIVE_HSIZE, scalarSpace);
    attNb.write(PredType::NATIVE_HSIZE, &nb);

    Attribute attDt = file.createAttribute("dt", PredType::NATIVE_DOUBLE, scalarSpace);
    attDt.write(PredType::NATIVE_DOUBLE, &dt);

    hsize_t nf = 0;
    Attribute attNf = file.createAttribute("num_frames", PredType::NATIVE_HSIZE, scalarSpace);
    attNf.write(PredType::NATIVE_HSIZE, &nf);

    {
        hsize_t initDims[1] = {nb};
        DataSpace initSpace(1, initDims);
        file.createDataSet("initial_data", GetInitType(), initSpace);
    }

    {
        hsize_t dims[2]    = {0, nb};
        hsize_t maxDims[2] = {H5S_UNLIMITED, nb};
        DataSpace tracksSpace(2, dims, maxDims);

        DSetCreatPropList prop;
        hsize_t chunkDims[2] = {1, nb}; 
        prop.setChunk(2, chunkDims);
        prop.setDeflate(4);

        file.createDataSet("tracks", GetPosType(), tracksSpace, prop);
    }

    file.flush(H5F_SCOPE_GLOBAL);
    return file;
}

void WriteSimulationFrame(H5::H5File& file, const std::vector<Object>& objs, std::size_t frameIndex)
{
    if (objs.empty()) return;
    hsize_t numBodies = static_cast<hsize_t>(objs.size());

    if (frameIndex == 0) {
        DataSet dset = file.openDataSet("initial_data");
        std::vector<InitRecord> buffer(objs.size());
        for (size_t i = 0; i < objs.size(); ++i) {
            buffer[i].position = objs[i].position;
            buffer[i].mass     = objs[i].mass;
            buffer[i].radius   = objs[i].radius;
        }
        dset.write(buffer.data(), GetInitType());
    }

    DataSet dset = file.openDataSet("tracks");
    
    hsize_t currentDims[2] = { static_cast<hsize_t>(frameIndex) + 1, numBodies };
    dset.extend(currentDims);

    DataSpace fileSpace = dset.getSpace();
    hsize_t offset[2] = { static_cast<hsize_t>(frameIndex), 0 };
    hsize_t count[2]  = { 1, numBodies };
    fileSpace.selectHyperslab(H5S_SELECT_SET, count, offset);

    std::vector<glm::dvec3> positions(objs.size());
    for (size_t i = 0; i < objs.size(); ++i) {
        positions[i] = objs[i].position;
    }

    hsize_t memDims[2] = { 1, numBodies };
    DataSpace memSpace(2, memDims);

    dset.write(positions.data(), GetPosType(), memSpace, fileSpace);
}

void CloseSimulationFile(H5::H5File& file, std::size_t totalFrames)
{
    try {
        Attribute att = file.openAttribute("num_frames");
        hsize_t nf = static_cast<hsize_t>(totalFrames);
        att.write(PredType::NATIVE_HSIZE, &nf);
        
        file.flush(H5F_SCOPE_GLOBAL);
        file.close();
    } catch (...) {
        std::cerr << "Error closing simulation file.\n";
    }
}
