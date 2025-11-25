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
        double getMass(); 
        std::size_t getBodiesAmount();
        glm::dvec3 getVel();
        glm::dvec3 getCenter(); 
        void transPointToSystem(std::vector<Object>& objects); // перенос точки отсчёта в центр масс системы 

    private: 
        // поля 
        std::size_t bodiesAmount; // количество объектов 
        glm::dvec3 systemVel; // вектор скорости системы 
        glm::dvec3 massCenter; // координаты центра масс
        double totalMass; // общая масса системы 

        //приватные методы 
        void calcParams(std::vector<Object>& objects); // подсчёт параметров системы 


}; 

#endif 