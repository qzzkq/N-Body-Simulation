#pragma once
#include <vector>
#include <string>
#include <memory>
#include <glm/glm.hpp>
#include "object.hpp"

struct GenParams {
    int count;
    double centralMass;
    double baseMass;
    float minRadius;
    float maxRadius;
    float spread;
    unsigned seed;
};

class IScenario {
public:
    virtual ~IScenario() = default;
    virtual std::string getName() const = 0;
    virtual void generate(std::vector<Object>& out, const GenParams& params) = 0;
};

class ScenarioManager {
public:
    void registerScenario(std::unique_ptr<IScenario> scenario);
    std::vector<std::string> getNames() const;
    void runScenario(size_t index, std::vector<Object>& out, const GenParams& params);
    bool isValidIndex(size_t index) const;

private:
    std::vector<std::unique_ptr<IScenario>> scenarios_;
};

std::unique_ptr<ScenarioManager> CreateDefaultManager();