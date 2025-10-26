#ifndef DATA_HPP
#define DATA_HPP
#pragma once

#include <vector>
#include <string>
#include <glm/glm.hpp>

// Формат записи: dvec3 position, dvec3 velocity, double mass, double radius
// (в файле хранится подряд 64 байта: 3*8 + 3*8 + 8 + 8)
struct Particle {
    glm::dvec3 position;
    glm::dvec3 velocity;
    double     mass;
    double     radius;
};

// Перезаписывает файл и пишет весь массив как 1D датасет с именем dsetName
void Writer(const std::string& fileName,
            const std::string& dsetName,
            const std::vector<Particle>& particles);

// Читает весь датасет целиком в вектор
std::vector<Particle> Reader(const std::string& fileName,
                             const std::string& dsetName);


#endif // DATA_HPP