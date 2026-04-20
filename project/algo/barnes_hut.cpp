// barnes_hut_step.cpp
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include "object.hpp"
#include "barnes_hut.hpp"
#include "physics.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif
struct OctreeNode {
    double mass;
    double cx, cy, cz;
    double x0, y0, z0, size;
    OctreeNode* children[8];
    Object* body;
    int leafCount;
    double leafMass;
    double leafMx, leafMy, leafMz;

    OctreeNode(double x, double y, double z, double s)
        : mass(0.0), cx(0.0), cy(0.0), cz(0.0),
          x0(x), y0(y), z0(z), size(s), body(nullptr),
          leafCount(0), leafMass(0.0), leafMx(0.0), leafMy(0.0), leafMz(0.0)
    {
        for (int i = 0; i < 8; ++i) children[i] = nullptr;
    }
};

static constexpr double MIN_NODE_SIZE_AU = 1.0e-8; // защита от бесконечного дробления при совпадающих координатах

static bool hasChildren(const OctreeNode* node) {
    for (int i = 0; i < 8; ++i) {
        if (node->children[i] != nullptr) {
            return true;
        }
    }
    return false;
}

static void clearLeafBodies(OctreeNode* node) {
    node->body = nullptr;
    node->leafCount = 0;
    node->leafMass = 0.0;
    node->leafMx = 0.0;
    node->leafMy = 0.0;
    node->leafMz = 0.0;
}

static void addLeafBody(OctreeNode* node, Object* body) {
    const auto p = body->GetPos();
    node->leafMass += body->mass;
    node->leafMx += p.x * body->mass;
    node->leafMy += p.y * body->mass;
    node->leafMz += p.z * body->mass;
    ++node->leafCount;
    node->body = (node->leafCount == 1) ? body : nullptr;
}

static int getOctant(const OctreeNode* node, double x, double y, double z) {
    double mx = node->x0 + node->size * 0.5;
    double my = node->y0 + node->size * 0.5;
    double mz = node->z0 + node->size * 0.5;
    int idx = 0;
    if (x >= mx) idx |= 1;
    if (y >= my) idx |= 2;
    if (z >= mz) idx |= 4;
    return idx;
}

static OctreeNode* createChildNode(OctreeNode* node,
                                   int childIdx,
                                   std::vector<OctreeNode>& pool,
                                   std::size_t poolMaxNodes) {
    if (node->children[childIdx] != nullptr) {
        return node->children[childIdx];
    }
    if (pool.size() >= poolMaxNodes) {
        return nullptr;
    }

    const double half = node->size * 0.5;
    const double nx = node->x0 + ((childIdx & 1) ? half : 0.0);
    const double ny = node->y0 + ((childIdx & 2) ? half : 0.0);
    const double nz = node->z0 + ((childIdx & 4) ? half : 0.0);
    pool.emplace_back(nx, ny, nz, half);
    node->children[childIdx] = &pool.back();
    return node->children[childIdx];
}

static bool insertBody(OctreeNode* node,
                       Object* body,
                       std::vector<OctreeNode>& pool,
                       std::size_t poolMaxNodes) {
    if (node->size <= MIN_NODE_SIZE_AU) {
        addLeafBody(node, body);
        return true;
    }

    auto p = body->GetPos();
    if (!hasChildren(node)) {
        if (node->leafCount == 0) {
            addLeafBody(node, body);
            return true;
        }

        if (node->leafCount > 1) {
            addLeafBody(node, body);
            return true;
        }

        Object* old = node->body;
        clearLeafBodies(node);

        const auto oldPos = old->GetPos();
        const int oldOct = getOctant(node, oldPos.x, oldPos.y, oldPos.z);
        OctreeNode* oldChild = createChildNode(node, oldOct, pool, poolMaxNodes);
        if (!oldChild || !insertBody(oldChild, old, pool, poolMaxNodes)) {
            for (auto& child : node->children) {
                child = nullptr;
            }
            clearLeafBodies(node);
            addLeafBody(node, old);
            addLeafBody(node, body);
            return true;
        }
    }

    const int idx = getOctant(node, p.x, p.y, p.z);
    OctreeNode* child = createChildNode(node, idx, pool, poolMaxNodes);
    if (!child) {
        return false;
    }
    return insertBody(child, body, pool, poolMaxNodes);
}

static void computeMass(OctreeNode* node) {
    if (!hasChildren(node)) {
        if (node->leafCount > 0 && node->leafMass > 0.0) {
            node->mass = node->leafMass;
            node->cx = node->leafMx / node->leafMass;
            node->cy = node->leafMy / node->leafMass;
            node->cz = node->leafMz / node->leafMass;
        } else {
            node->mass = 0.0;
        }
        return;
    }
    double m = 0.0, mx = 0.0, my = 0.0, mz = 0.0;
    for (int i = 0; i < 8; ++i) {
        OctreeNode* c = node->children[i];
        if (!c) continue;
        computeMass(c);
        if (c->mass <= 0.0) continue;
        m  += c->mass;
        mx += c->cx * c->mass;
        my += c->cy * c->mass;
        mz += c->cz * c->mass;
    }
    node->mass = m;
    if (m > 0.0) {
        node->cx = mx / m;
        node->cy = my / m;
        node->cz = mz / m;
    }
}

static void accumulateAccel(const Object& obj, const OctreeNode* node,
                            double px, double py, double pz,
                            double& ax, double& ay, double& az,
                            double theta,
                            double softeningSq)
{
    if (node->mass <= 0.0) return;

    if (!hasChildren(node)) {
        if (node->leafCount == 1 && node->body == &obj) return;
        double dx = node->cx - px;
        double dy = node->cy - py;
        double dz = node->cz - pz;
        double r2 = dx*dx + dy*dy + dz*dz + softeningSq;
        if (r2 <= 0.0) return;
        double invR = 1.0 / std::sqrt(r2);
        double invR3 = invR * invR * invR;
        double f = physics::G * node->mass * invR3;
        ax += f * dx;
        ay += f * dy;
        az += f * dz;
        return;
    }

    double dx = node->cx - px;
    double dy = node->cy - py;
    double dz = node->cz - pz;
    double r2 = dx*dx + dy*dy + dz*dz + softeningSq;
    if (r2 <= 0.0) return;
    double r = std::sqrt(r2);

    if ((node->size / r) < theta) {
        double invR = 1.0 / r;
        double invR3 = invR * invR * invR;
        double f = physics::G * node->mass * invR3;
        ax += f * dx;
        ay += f * dy;
        az += f * dz;
    } else {
        for (int i = 0; i < 8; ++i) {
            if (node->children[i]) {
                accumulateAccel(obj, node->children[i], px, py, pz, ax, ay, az, theta, softeningSq);
            }
        }
    }
}

void simulationStepBarnesHutCPU(std::vector<Object>& objs, double dt, bool pause, bool forceSync) {
    (void)forceSync;
    if (pause || objs.empty()) return;
    if (dt <= 0.0) return;

    const std::size_t n = objs.size();
    const double theta = physics::getBarnesHutTheta();
    const double soft = physics::getSofteningAU();
    const double softSq = soft * soft;
    const std::size_t poolMaxNodes = std::max<std::size_t>(100000, n * 256ull);
    static std::vector<OctreeNode> pool;
    if (pool.capacity() < poolMaxNodes) {
        pool.reserve(poolMaxNodes);
    }
    static std::vector<Object> temp;
    static std::vector<double> masses;
    static std::vector<glm::dvec3> x, v;
    static std::vector<glm::dvec3> k1x, k2x, k3x, k4x;
    static std::vector<glm::dvec3> k1v, k2v, k3v, k4v;
    static std::vector<glm::dvec3> tmpX, tmpV;
    temp = objs;
    masses.resize(n);
    x.resize(n); v.resize(n);
    k1x.resize(n); k2x.resize(n); k3x.resize(n); k4x.resize(n);
    k1v.resize(n); k2v.resize(n); k3v.resize(n); k4v.resize(n);
    tmpX.resize(n); tmpV.resize(n);

    for (std::size_t i = 0; i < n; ++i) {
        masses[i] = objs[i].mass;
        x[i] = objs[i].position;
        v[i] = objs[i].velocity;
    }

    auto computeAccelerationsDirect = [&](const std::vector<glm::dvec3>& positions,
                                          std::vector<glm::dvec3>& outAcc,
                                          double* outMaxAcc = nullptr) {
        outAcc.assign(n, glm::dvec3(0.0));
        double maxAccGlobal = 0.0;
#ifdef _OPENMP
        #pragma omp parallel for reduction(max:maxAccGlobal) schedule(dynamic, 64)
#endif
        for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(n); ++i) {
            glm::dvec3 ai(0.0);
            for (std::size_t j = 0; j < n; ++j) {
                if (static_cast<std::size_t>(i) == j) continue;
                const glm::dvec3 delta = positions[j] - positions[static_cast<std::size_t>(i)];
                const double rawDistSq = delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
                const double distSq = rawDistSq + softSq;
                if (rawDistSq <= 0.0 && softSq <= 0.0) {
                    continue;
                }
                const double invDist = 1.0 / std::sqrt(distSq);
                const double invDist3 = invDist * invDist * invDist;
                ai += delta * (physics::G * invDist3 * masses[j]);
            }
            outAcc[static_cast<std::size_t>(i)] = ai;
            maxAccGlobal = std::max(maxAccGlobal, glm::length(ai));
        }
        if (outMaxAcc != nullptr) {
            *outMaxAcc = maxAccGlobal;
        }
    };

    auto computeAccelerationsBH = [&](const std::vector<glm::dvec3>& positions,
                                      std::vector<glm::dvec3>& outAcc,
                                      double* outMaxAcc = nullptr) {
        outAcc.assign(n, glm::dvec3(0.0));
        for (std::size_t i = 0; i < n; ++i) {
            temp[i].position = positions[i];
            temp[i].mass = masses[i];
        }

        auto p0 = temp[0].GetPos();
        double minX = p0.x, maxX = p0.x;
        double minY = p0.y, maxY = p0.y;
        double minZ = p0.z, maxZ = p0.z;

        for (const auto& o : temp) {
            auto p = o.GetPos();
            if (p.x < minX) minX = p.x;
            if (p.x > maxX) maxX = p.x;
            if (p.y < minY) minY = p.y;
            if (p.y > maxY) maxY = p.y;
            if (p.z < minZ) minZ = p.z;
            if (p.z > maxZ) maxZ = p.z;
        }

        double spanX = maxX - minX;
        double spanY = maxY - minY;
        double spanZ = maxZ - minZ;
        double side = std::max({spanX, spanY, spanZ});
        if (side <= 0.0) side = 1.0;

        double cx = 0.5 * (minX + maxX);
        double cy = 0.5 * (minY + maxY);
        double cz = 0.5 * (minZ + maxZ);

        pool.clear();
        pool.emplace_back(cx - side * 0.5, cy - side * 0.5, cz - side * 0.5, side);
        OctreeNode* root = &pool.back();
        bool buildOk = true;
        for (auto& t : temp) {
            if (!insertBody(root, &t, pool, poolMaxNodes)) {
                buildOk = false;
                break;
            }
        }
        if (!buildOk) {
            computeAccelerationsDirect(positions, outAcc, outMaxAcc);
            return;
        }
        computeMass(root);

        double maxA = 0.0;
#ifdef _OPENMP
        #pragma omp parallel for reduction(max:maxA) schedule(dynamic, 128)
#endif
        for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(n); ++i) {
            auto p = temp[static_cast<std::size_t>(i)].GetPos();
            double ax = 0.0, ay = 0.0, az = 0.0;
            accumulateAccel(temp[static_cast<std::size_t>(i)], root, p.x, p.y, p.z, ax, ay, az, theta, softSq);
            outAcc[static_cast<std::size_t>(i)] = glm::dvec3(ax, ay, az);
            maxA = std::max(maxA, glm::length(outAcc[static_cast<std::size_t>(i)]));
        }
        if (outMaxAcc) *outMaxAcc = maxA;
    };

    double remaining = dt;
    constexpr int MAX_SUBSTEPS = 1024;
    const double finishEps = std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(dt));
    int substeps = 0;
    while (remaining > 0.0 && substeps < MAX_SUBSTEPS) {
        double maxAcc = 0.0;
        if (substeps == 0) {
            computeAccelerationsBH(x, k1v, &maxAcc); // a(t)
        } else {
            maxAcc = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                maxAcc = std::max(maxAcc, glm::length(k1v[i]));
            }
        }

        double suggested = remaining;
        if (maxAcc > 0.0) {
            suggested = std::min(suggested, 0.02 / std::sqrt(maxAcc));
        }

        const int stepsLeft = MAX_SUBSTEPS - substeps;
        const double minStep = remaining / static_cast<double>(stepsLeft);
        double h = std::clamp(suggested, minStep, remaining);

        // Velocity Verlet
        // 1. x(t+h) = x(t) + v(t)*h + 0.5*a(t)*h^2
        for (std::size_t i = 0; i < n; ++i) {
            tmpX[i] = x[i] + v[i] * h + 0.5 * k1v[i] * h * h;
        }

        // 2. a(t+h)
        computeAccelerationsBH(tmpX, k2v, nullptr);

        // 3. v(t+h) = v(t) + 0.5*(a(t) + a(t+h))*h
        for (std::size_t i = 0; i < n; ++i) {
            v[i] += 0.5 * h * (k1v[i] + k2v[i]);
            x[i] = tmpX[i];
            k1v[i] = k2v[i];
        }

        remaining = std::max(0.0, remaining - h);
        if (remaining <= finishEps) {
            remaining = 0.0;
        }
        ++substeps;
    }

    if (substeps == 0) {
        computeAccelerationsBH(x, k1v, nullptr);
    }
    
    for (std::size_t i = 0; i < n; ++i) {
        objs[i].position = x[i];
        objs[i].velocity = v[i];
        objs[i].acceleration = k1v[i]; // k1v contains the final acceleration
    }
}
