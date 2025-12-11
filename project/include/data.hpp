#ifndef DATA_HPP
#define DATA_HPP
#pragma once

#include <vector>
#include <string>
#include <glm/glm.hpp>

class Object;

// Вперё объявляем H5::H5File, чтобы не тащить H5Cpp.h в заголовок
namespace H5 {
    class H5File;
}

struct Particle {
    glm::dvec3 position;
    glm::dvec3 velocity;
    double     mass;
    double     radius;
};

void Writer(const std::string& fileName,
            const std::string& dsetName,
            const std::vector<Particle>& particles);

std::vector<Particle> Reader(const std::string& fileName,
                             const std::string& dsetName);

std::vector<std::string> ListH5Files(const std::string& dir = "data");

bool LoadObjectsFromFile(const std::string& filePath,
                         const std::string& dsetName,
                         std::vector<Object>& outObjs);

H5::H5File OpenFramesFile(const std::string& fileName,
                          std::size_t numBodies);

void WriteFrame(H5::H5File& file,
                const std::vector<Object>& objs,
                double t,
                std::size_t frameIndex,
                const std::string& prefix = "frame");

void FinalizeFramesFile(H5::H5File& file,
                        std::size_t numFrames);

#endif // DATA_HPP
