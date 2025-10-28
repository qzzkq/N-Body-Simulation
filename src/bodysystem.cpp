#include <glm/glm.hpp>
#include <vector>
#include "object.hpp" 
#include "bodysystem.hpp"

BodySystem::BodySystem(std::vector<Object>& objects) : bodies_(objects) {

    this->bodies_ = objects; 
    this->bodiesAmount = objects.size();
    double tMass = 0.;
    glm::dvec3 impulse = glm::dvec3(0., 0., 0.);  

    if (this->bodiesAmount == 0) {
        this->systemVel = impulse; 
        this->totalMass = 0; 
    }
    
    for (size_t i = 0; i < this->bodiesAmount; ++i) {
        tMass += this->bodies_.at(i).mass;
        impulse += this->bodies_.at(i).mass * this->bodies_.at(i).velocity; 
    }

    this->systemVel = impulse / tMass; 
    this->totalMass = tMass; 

}

glm::dvec3 BodySystem::getVel() {

    glm::dvec3 impulse = glm::dvec3(0., 0., 0.); 

    if (this->bodiesAmount == 0) {
        return  glm::dvec3(0., 0., 0.); 
    }
    
    for (size_t i = 0; i < this->bodiesAmount; ++i) {
        impulse += this->bodies_.at(i).mass * this->bodies_.at(i).velocity; 
    }

    this->systemVel = impulse / this->totalMass; 

    return this->systemVel;
} 
double BodySystem::getMass() {
    return this->totalMass; 
} 
std::size_t BodySystem::getBodiesAmount() {
    return this->bodiesAmount; 
}

