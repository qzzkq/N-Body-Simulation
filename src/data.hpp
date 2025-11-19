#ifndef DATA_HPP
#define DATA_HPP
#pragma once

#include <vector>
#include <string>
#include <glm/glm.hpp>
#include <H5Cpp.h>

class Object;

struct Particle {
    glm::dvec3 position;
    glm::dvec3 velocity;
    double     mass;
    double     radius;
};

// Запись одного датасета частиц в HDF5 (как снимок системы)
void Writer(const std::string& fileName,
            const std::string& dsetName,
            const std::vector<Particle>& particles);

// Чтение датасета целиком
std::vector<Particle> Reader(const std::string& fileName,
                             const std::string& dsetName);

// Список .h5 файлов в каталоге (по умолчанию "data")
std::vector<std::string> ListH5Files(const std::string& dir = "data");

// Грузим объекты из файла (через Reader)
bool LoadObjectsFromFile(const std::string& filePath,
                         const std::string& dsetName,
                         std::vector<Object>& outObjs);

// Открыть файл для серии кадров (затирает старый файл)
H5::H5File OpenFramesFile(const std::string& fileName,
                          std::size_t numBodies);

// Записать один кадр: frame_000000, frame_000001, и тд
void WriteFrame(H5::H5File& file,
                const std::vector<Object>& objs,
                double t,
                std::size_t frameIndex,
                const std::string& prefix = "frame");

// В конце симуляции обновить атрибут num_frames
void FinalizeFramesFile(H5::H5File& file,
                        std::size_t numFrames);

#endif // DATA_HPP