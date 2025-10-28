#pragma once 
#ifndef SYSTEMINC
#define SYSTEMINC

#include <glm/glm.hpp>
#include <vector>
#include "object.hpp" 

class BodySystem {

    public: 
        BodySystem(std::vector<Object>& objects);
        // получение информации о системе  
        glm::dvec3 getVel(); 
        double getMass(); 
        std::size_t getBodiesAmount(); 

    private: 
        std::vector<Object>& bodies_; // ссылка на имеющиеся объекты в системе 
        std::size_t bodiesAmount; // количество объектов 
        glm::dvec3 systemVel; // вектор скорости системы 
        double totalMass; // общая масса системы 

}; 

#endif 