#include "physics.hpp"

namespace physics {

void colorFromMass(std::vector<Object>& objs) {
    if (objs.empty()) return;

    double minMass = std::numeric_limits<double>::max();
    double maxMass = -1.0;

    // определение минимальной массы 
    for (const auto& obj : objs) {
        if (obj.mass > maxMass) maxMass = obj.mass;
        if (obj.mass < minMass) minMass = obj.mass;
    }

    // Защита от нулевой массы и log(0)
    minMass = std::max(minMass, 1e-30);
    maxMass = std::max(maxMass, 1e-30);

    // 2. Палитра цветов для планет 
    const glm::vec3 colorCold = glm::vec3(0.8f, 0.4f, 0.0f); // Темно-оранжевый (лёгкие)
    const glm::vec3 colorMid  = glm::vec3(1.0f, 0.9f, 0.1f); // Ярко-желтый (средние по массе)
    const glm::vec3 colorHot  = glm::vec3(1.0f, 1.0f, 1.0f); // Белый (самые тяжёлые)

    // Используем логарифмическую шкалу, так как массы планет отличаются в миллиарды раз
    double logMin = std::log10(minMass);
    double logMax = std::log10(maxMass);
    double range = logMax - logMin;

    // присваиваем цвета объектам
    for (auto& obj : objs) {
        double m = std::max(obj.mass, 1e-30);
        double val = std::log10(m);
        float t = 0.0f;
        if (range > 1e-9) {
            t = static_cast<float>((val - logMin) / range);
        }
        t = std::clamp(t, 0.0f, 1.0f);
        glm::vec3 finalColor;
        if (t < 0.6f) {
            float localT = t / 0.6f; 
            finalColor = glm::mix(colorCold, colorMid, localT);
        } 
        else {
            float localT = (t - 0.6f) / 0.4f;
            finalColor = glm::mix(colorMid, colorHot, localT);
        }
        obj.color = glm::vec4(finalColor, 1.0f);
    }
}

}