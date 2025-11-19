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
        std::cerr << "HDF5 Reader error: " << e.getDetailMsg()
                  << " (file=" << fileName
                  << ", dset=" << dsetName << ")\n";
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
                 1410.0f,           // density
                 std::nullopt);     // radius пока не задаём явно

        o.Initalizing = false;

        // Если в файле есть валидный радиус — используем его.
        if (std::isfinite(p.radius) && p.radius > 0.0) {
            o.radius = static_cast<float>(p.radius);
        } else {
            // Иначе пробуем вычислить из массы и плотности.
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

    return file;
}

void WriteFrame(H5::H5File& file,
                const std::vector<Object>& objs,
                double t,
                std::size_t frameIndex,
                const std::string& prefix)
{
    if (objs.empty()) return;
    std::vector<Particle> parts;
    parts.reserve(objs.size());
    for (const auto& o : objs) {
        Particle p;
        p.position = o.position;
        p.velocity = o.velocity;
        p.mass     = o.mass;
        p.radius   = o.radius;
        parts.push_back(p);
    }

    H5::CompType& type = GetParticleType();

    hsize_t dims[1] = { static_cast<hsize_t>(parts.size()) };
    H5::DataSpace space(1, dims);

    char name[64];
    std::snprintf(name, sizeof(name), "%s_%06zu",
                  prefix.c_str(), frameIndex);

    H5::DataSet dset = file.createDataSet(name, type, space);

    dset.write(parts.data(), type);

    // Атрибут времени кадра t
    H5::DataSpace scalar(H5S_SCALAR);
    H5::Attribute tAttr = dset.createAttribute(
        "t",
        H5::PredType::NATIVE_DOUBLE,
        scalar
    );
    tAttr.write(H5::PredType::NATIVE_DOUBLE, &t);

    // Сразу выталкиваем данные на диск
    file.flush(H5F_SCOPE_GLOBAL);
}

void FinalizeFramesFile(H5::H5File& file,
                        std::size_t numFrames)
{
    // Перезаписываем num_frames
    H5::Attribute aNf = file.openAttribute("num_frames");
    hsize_t nf = static_cast<hsize_t>(numFrames);
    aNf.write(H5::PredType::NATIVE_HSIZE, &nf);

    file.flush(H5F_SCOPE_GLOBAL);
}