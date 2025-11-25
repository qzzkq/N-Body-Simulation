#ifndef OBJECT_HPP
#define OBJECT_HPP

#include <glm/glm.hpp>
#include <optional>
#include <vector>

class Object {
public:
    // constructor 
    Object(glm::dvec3 initPosition, 
        glm::dvec3 initVelocity, 
        double mass, 
        glm::dvec3 acc);


    //getters 
    glm::dvec3 getPos() const;
    glm::dvec3 getVel() const;
    glm::dvec3 getLastPos() const;
    glm::dvec3 getAcc() const;
    glm::vec4 getColor() const;
    double getMass() const; 

    //setters
    void setMass(double newMass);
    void setPos(glm::dvec3); 
    void setVel(glm::dvec3); 
    void setLastPos(glm::dvec3);
    void setAcc(glm::dvec3); 
    void setColor(glm::vec4); 

    bool Initalizing = false;
    bool Launched = false;

private:
    glm::dvec3 position; // current position 
    glm::dvec3 velocity; // object velocity  
    glm::dvec3 lastPos; // last step position
    glm::dvec3 acceleration; // acceleration 
    glm::vec4 color; // object color 
    double mass; // object mass 
};

#endif // OBJECT_HPP