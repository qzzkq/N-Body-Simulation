#include <glm/glm.hpp>
#include <vector>
#include "object.hpp" 
#include "bodysystem.hpp"

BodySystem::BodySystem(std::vector<Object>& objects) {
    calcParams(objects); 
}

double BodySystem::getMass() {
    return this->totalMass; 
} 

std::size_t BodySystem::getBodiesAmount() {
    return this->bodiesAmount; 
}

glm::dvec3 BodySystem::getVel() {
    return this->systemVel; 
}

glm::dvec3 BodySystem::getCenter() {
    return this->massCenter; 
} 

void BodySystem::calcParams(std::vector<Object>& objects) {

    glm::dvec3 impulse = glm::dvec3(0., 0., 0.); 
    glm::dvec3 center = glm::dvec3(0., 0., 0.);
    glm::dvec3 massMoment = glm::dvec3(0., 0., 0.);
    double tMass = 0.; 
    this->bodiesAmount = objects.size(); 

    if (this->bodiesAmount == 0) {
        this->systemVel = impulse; 
        this->totalMass = 0;
        this->massCenter = center; 
        return; 
    }
    
    for (size_t i = 0; i < this->bodiesAmount; ++i) {
        tMass += objects.at(i).GetMass(); // вычисляем массу системы 
        impulse += objects.at(i).GetMass() * objects.at(i).GetVel(); // вычисляем импульс 
        massMoment += objects.at(i).GetMass() * objects.at(i).GetPos(); // сумма произведений массы на координаты 

    }

    // считаем центр масс 
    center = massMoment / tMass; 

    // записываем параметры 
    this->massCenter = center; 
    this->totalMass = tMass; 
    this->systemVel = impulse / this->totalMass; 

}

void BodySystem::transPointToSystem(std::vector<Object>& objects) {

    for (size_t i = 0; i < this->bodiesAmount; ++i) {
        objects.at(i).GetVel() -= this->systemVel; 
    }

}

