#ifndef DATA_HPP
#define DATA_HPP
#pragma once

#include <vector>
#include <string>
#include <glm/glm.hpp>
#include "H5Cpp.h"
#include "object.hpp"
#include "graphic_state.hpp"

struct Particle {
    glm::dvec3 position;
    glm::dvec3 velocity;
    double     mass;
    double     radius;
    glm::vec4  color{1.0f, 1.0f, 1.0f, 1.0f};
};

H5::H5File CreateSimulationFile(const std::string& fileName, std::size_t numBodies, double dt);

void WriteSimulationFrame(H5::H5File& file, const std::vector<Object>& objs,
                          const std::vector<GraphicState>& graphics, std::size_t frameIndex);

void CloseSimulationFile(H5::H5File& file, std::size_t totalFrames);

void Writer(const std::string& fileName,
            const std::string& dsetName,
            const std::vector<Particle>& particles);

std::vector<Particle> Reader(const std::string& fileName,
                             const std::string& dsetName,
                             bool* outFileHadColorMember = nullptr);

std::vector<std::string> ListH5Files(const std::string& dir = "data");

bool LoadObjectsFromFile(const std::string& filePath,
                         const std::string& dsetName,
                         std::vector<Object>& outObjs,
                         std::vector<GraphicState>* outGraphics = nullptr,
                         bool* outHadColorInFile = nullptr);

bool LoadSystemFromTextFile(const std::string& filePath,
                            std::vector<Object>& outObjs,
                            std::vector<GraphicState>* outGraphics = nullptr);

bool SaveSystemToTextFile(const std::string& filePath,
                          const std::vector<Object>& objs);

#endif
