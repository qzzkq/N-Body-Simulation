#include "brut_force.hpp"
#include <cmath>
#include <algorithm>
#include <object.hpp>
#include <physics.hpp>
float G = physics::G;
void simulationStepBrutForceCPU(std::vector<Object>& objs, float dt, bool pause, int iterations) {
    for (int iter = 0; iter < iterations; ++iter) {
        for (size_t i = 0; i < objs.size(); ++i) {
            Object& obj = objs[i];

            for (size_t j = i + 1; j < objs.size(); ++j) {
                Object& obj2 = objs[j];

                glm::dvec3 delta = obj2.GetPos() - obj.GetPos();
                double distance = std::sqrt(delta.x * delta.x + delta.y * delta.y + delta.z * delta.z);
                if (distance <= 0.0) {
                    continue;
                }

                glm::dvec3 dir = delta / distance;
                double combinedRadius = obj.radius + obj2.radius;

                double effectiveDistance = std::max(distance, combinedRadius);
                double dist_m = effectiveDistance * 1000.0;        // м
                double F = (G * obj.mass * obj2.mass) / (dist_m * dist_m);              // Н
                float acc1_kmps2 = static_cast<float>((F / obj.mass)  / 1000.0);        // км/с²
                float acc2_kmps2 = static_cast<float>((F / obj2.mass) / 1000.0);        // км/с²
                glm::vec3 accObj  = glm::vec3(dir * static_cast<double>(acc1_kmps2));
                glm::vec3 accObj2 = glm::vec3(-dir * static_cast<double>(acc2_kmps2));

                if (!pause) {
                    obj.accelerate(accObj.x, accObj.y, accObj.z, dt);
                    obj2.accelerate(accObj2.x, accObj2.y, accObj2.z, dt);
                }
                /*
                if (distance < combinedRadius) {
                    glm::vec3 normal = glm::vec3(dir);
                    glm::vec3 relativeVelocity = glm::vec3(obj.velocity - obj2.velocity);
                    float relVelAlongNormal = glm::dot(relativeVelocity, normal);

                    if (relVelAlongNormal < 0.0f) {
                        double restitution = 0.8;
                        double invMass1 = 1.0 / obj.mass;
                        double invMass2 = 1.0 / obj2.mass;
                        double impulseScalar = -(1.0 + restitution) * static_cast<double>(relVelAlongNormal) / (invMass1 + invMass2);
                        glm::dvec3 impulse = glm::dvec3(normal) * impulseScalar;
                        obj.velocity += impulse * invMass1;
                        obj2.velocity -= impulse * invMass2;
                    }

                    double penetration = combinedRadius - distance;
                    if (penetration > 0.0) {
                        double invMass1 = 1.0 / obj.mass;
                        double invMass2 = 1.0 / obj2.mass;
                        double invMassSum = invMass1 + invMass2;
                        if (invMassSum > 0.0) {
                            double correctionScale = penetration / invMassSum;
                            glm::dvec3 correction = glm::dvec3(normal) * correctionScale;
                            obj.position -= correction * invMass1;
                            obj2.position += correction * invMass2;
                        }
                    }
                */
            }
            if (!pause) {
                obj.UpdatePos(dt);
            }
        }
    }
}