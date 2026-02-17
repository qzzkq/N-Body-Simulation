#ifndef SIMSTATEINC
#define SIMSTATEINC

typedef struct SimState {
    bool running = true;
    bool pause = false;
    float timeScale = 1.0f;
    float deltaTime = 0.0f;
    float lastFrame = 0.0f;
} SimState;

#endif 