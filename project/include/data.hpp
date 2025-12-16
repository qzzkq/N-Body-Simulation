#ifndef DATA_HPP
#define DATA_HPP
#pragma once

#include <vector>
#include <string>
#include <glm/glm.hpp>
#include "H5Cpp.h"
#include "object.hpp"

struct Particle {
    glm::dvec3 position;
    glm::dvec3 velocity;
    double     mass;
    double     radius;
};

H5::H5File CreateSimulationFile(const std::string& fileName, std::size_t numBodies, double dt);

void WriteSimulationFrame(H5::H5File& file, const std::vector<Object>& objs, std::size_t frameIndex);

void CloseSimulationFile(H5::H5File& file, std::size_t totalFrames);

void Writer(const std::string& fileName,
            const std::string& dsetName,
            const std::vector<Particle>& particles);

std::vector<Particle> Reader(const std::string& fileName,
                             const std::string& dsetName);

std::vector<std::string> ListH5Files(const std::string& dir = "data");

bool LoadObjectsFromFile(const std::string& filePath,
                         const std::string& dsetName,
                         std::vector<Object>& outObjs);

#endif