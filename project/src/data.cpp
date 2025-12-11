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

        o.UpdateVertices();
        outObjs.push_back(std::move(o));
    }

    return true;
}


// frame_000000: position + mass + radius
// frame_000001+: только position

struct Frame0Record {
    glm::dvec3 position;
    double     mass;
    double     radius;
};

static CompType MakeFrame0Type()
{
    hsize_t v3dims[1] = {3};
    ArrayType vec3Type(PredType::NATIVE_DOUBLE, 1, v3dims);

    CompType t(sizeof(Frame0Record));
    t.insertMember("position", HOFFSET(Frame0Record, position), vec3Type);
    t.insertMember("mass",     HOFFSET(Frame0Record, mass),     PredType::NATIVE_DOUBLE);
    t.insertMember("radius",   HOFFSET(Frame0Record, radius),   PredType::NATIVE_DOUBLE);
    return t;
}

static CompType& GetFrame0Type()
{
    static CompType t = MakeFrame0Type();
    return t;
}

H5::H5File OpenFramesFile(const std::string& fileName,
                          std::size_t numBodies)
{
    H5::H5File file(fileName, H5F_ACC_TRUNC);

    hsize_t dims[1] = {1};
    H5::DataSpace scalar(1, dims);

    {
        hsize_t nb = static_cast<hsize_t>(numBodies);
        H5::Attribute aNb = file.createAttribute(
            "num_bodies",
            H5::PredType::NATIVE_HSIZE,
            scalar
        );
        aNb.write(H5::PredType::NATIVE_HSIZE, &nb);
    }

    {
        hsize_t nf = 0;
        H5::Attribute aNf = file.createAttribute(
            "num_frames",
            H5::PredType::NATIVE_HSIZE,
            scalar
        );
        aNf.write(H5::PredType::NATIVE_HSIZE, &nf);
    }

    file.flush(H5F_SCOPE_GLOBAL);
    return file;
}

void WriteFrame(H5::H5File& file,
                const std::vector<Object>& objs,
                double t,
                std::size_t frameIndex,
                const std::string& prefix)
{
    if (objs.empty()) return;

    char name[64];
    std::snprintf(name, sizeof(name), "%s_%06zu",
                  prefix.c_str(), static_cast<std::size_t>(frameIndex));

    hsize_t dims[1] = { static_cast<hsize_t>(objs.size()) };
    H5::DataSpace space(1, dims);

    if (frameIndex == 0) {
        std::vector<Frame0Record> recs;
        recs.reserve(objs.size());
        for (const auto& o : objs) {
            Frame0Record r;
            r.position = o.position;
            r.mass     = o.mass;
            r.radius   = o.radius;
            recs.push_back(r);
        }

        H5::CompType& type = GetFrame0Type();
        H5::DataSet dset = file.createDataSet(name, type, space);
        dset.write(recs.data(), type);

        H5::DataSpace scalar(H5S_SCALAR);
        H5::Attribute tAttr = dset.createAttribute(
            "time",
            H5::PredType::NATIVE_DOUBLE,
            scalar
        );
        tAttr.write(H5::PredType::NATIVE_DOUBLE, &t);
    } else {
        std::vector<glm::dvec3> pos;
        pos.reserve(objs.size());
        for (const auto& o : objs) {
            pos.push_back(o.position);
        }

        hsize_t v3dims[1] = {3};
        H5::ArrayType vec3Type(PredType::NATIVE_DOUBLE, 1, v3dims);

        H5::DataSet dset = file.createDataSet(name, vec3Type, space);
        dset.write(pos.data(), vec3Type);

        H5::DataSpace scalar(H5S_SCALAR);
        H5::Attribute tAttr = dset.createAttribute(
            "time",
            H5::PredType::NATIVE_DOUBLE,
            scalar
        );
        tAttr.write(H5::PredType::NATIVE_DOUBLE, &t);
    }

    file.flush(H5F_SCOPE_GLOBAL);
}

void FinalizeFramesFile(H5::H5File& file,
                        std::size_t numFrames)
{
    try {
        H5::Attribute aNf = file.openAttribute("num_frames");
        hsize_t nf = static_cast<hsize_t>(numFrames);
        aNf.write(H5::PredType::NATIVE_HSIZE, &nf);
        file.flush(H5F_SCOPE_GLOBAL);
    } catch (const H5::Exception& e) {
        std::cerr << "FinalizeFramesFile error: " << e.getDetailMsg() << "\n";
    }
}
