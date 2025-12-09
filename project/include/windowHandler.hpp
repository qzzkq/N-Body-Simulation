#ifndef WHANDLER_HPP
#define WHANDLER_HPP
#pragma once
#include <GL/glew.h>
#include <GLFW/glfw3.h>

class WindowHandler {

    public: 

        WindowHandler(); 

        double getTime(); // получить время в симуляции

        int createWindow(); // создать окно 
        
        void setTitle(char* str); // вывести информацию в окно 

        void cleanResources(); // очистка ресурсов 

        void show(); // показать кадр 

        GLFWwindow* getWindowPointer(); // получить указатель на окно

        int closeFlag(); // флаг закрытия окна 

        void cleanWindow(); // очистить буферы (по факту окно)

    private: 

        GLFWwindow* window; 

};

#endif