#include "object.hpp"
#include <glm/gtc/constants.hpp>
#include <cmath>


Object::Object(glm::dvec3 initPosition, 
        glm::dvec3 initVelocity, 
        double mass, 
        glm::dvec3 acc)
{
    this->position = initPosition;
    this->velocity = initVelocity;
    this->color = glm::vec4(1.0f, 0.0f, 0.0f, 1.0f);
    this->lastPos = this->position;
    this->mass = mass; 
    this->acceleration = glm::dvec3(0.0);
}

//getters 
glm::dvec3 Object::getPos() const {
    return this->position;
}
glm::dvec3 Object::getVel() const {
    return this->velocity; 
}
glm::dvec3 Object::getLastPos() const {
    return this->lastPos; 
}
glm::dvec3 Object::getAcc() const {
    return this->acceleration; 
}
glm::vec4 Object::getColor() const {
    return this->color; 
}
double Object::getMass() const {
    return this->mass; 
}

//setters
void Object::setMass(double newMass) {
    this->mass = newMass; 
}
void Object::setPos(glm::dvec3 newPos) {
    this->position = newPos; 
} 
void Object::setVel(glm::dvec3 newVel) {
    this->velocity = newVel; 
} 
void Object::setLastPos(glm::dvec3 newLastPos) {
    this->lastPos = newLastPos; 
}
void Object::setAcc(glm::dvec3 newAcc) {
    this->acceleration = newAcc; 
} 
void Object::setColor(glm::vec4 newColor) {
    this->color = newColor; 
}