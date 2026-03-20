// barnes_hut_step.cpp
#include <vector>
#include <cmath>
#include <algorithm>
#include "object.hpp"
#include "barnes_hut.hpp"
#include "physics.hpp"
struct OctreeNode {
    double mass;
    double cx, cy, cz;
    double x0, y0, z0, size;
    OctreeNode* children[8];
    Object* body;

    OctreeNode(double x, double y, double z, double s)
        : mass(0.0), cx(0.0), cy(0.0), cz(0.0),
          x0(x), y0(y), z0(z), size(s), body(nullptr)
    {
        for (int i = 0; i < 8; ++i) children[i] = nullptr;
    }
};

static constexpr double MIN_NODE_SIZE_AU = 1.0e-8; // защита от бесконечного дробления при совпадающих координатах

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

static void insertBody(OctreeNode* node,
                        Object* body,
                        std::vector<OctreeNode>& pool,
                        std::size_t poolMaxNodes) {
    // Если узел уже "слишком мал", объединяем тела в одном листе.
    // (Иначе при совпадающих координатах рекурсия будет дробить дерево бесконечно.)
    if (node->size <= MIN_NODE_SIZE_AU) {
        if (!node->body) node->body = body; // второе тело игнорируем
        return;
    }

    auto createNode = [&](double nx, double ny, double nz, double s) -> OctreeNode* {
        if (pool.size() >= poolMaxNodes) return nullptr;
        pool.emplace_back(nx, ny, nz, s);
        return &pool.back();
    };

    auto p = body->GetPos();
    if (!node->children[0] && !node->body) {
        node->body = body;
        return;
    }
    if (!node->children[0] && node->body) {
        Object* old = node->body;
        node->body = nullptr;
        double half = node->size * 0.5;
        bool anyChild = false;
        for (int i = 0; i < 8; ++i) {
            double nx = node->x0 + ((i & 1) ? half : 0.0);
            double ny = node->y0 + ((i & 2) ? half : 0.0);
            double nz = node->z0 + ((i & 4) ? half : 0.0);
            node->children[i] = createNode(nx, ny, nz, half);
            anyChild = anyChild || (node->children[i] != nullptr);
        }
        if (!anyChild) {
            // Пул исчерпан: не дробим дерево, оставляем старое тело в узле.
            node->body = old;
            return;
        }

        int o = getOctant(node, old->GetPos().x, old->GetPos().y, old->GetPos().z);
        if (node->children[o]) {
            node->children[o]->body = old;
        } else {
            // Не можем поместить старое тело в нужный октант из-за ограничения пула.
            node->body = old;
            for (int i = 0; i < 8; ++i) node->children[i] = nullptr;
            return;
        }
    }
    int idx = getOctant(node, p.x, p.y, p.z);
    OctreeNode* child = node->children[idx];
    if (!child) return;
    if (!child->children[0] && !child->body) {
        child->body = body;
    } else {
        insertBody(child, body, pool, poolMaxNodes);
    }
}

static void computeMass(OctreeNode* node) {
    if (!node->children[0]) {
        if (node->body) {
            node->mass = node->body->mass;
            auto p = node->body->GetPos();
            node->cx = p.x;
            node->cy = p.y;
            node->cz = p.z;
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

    if (!node->children[0]) {
        if (node->body == &obj) return;
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
    const std::size_t poolMaxNodes = std::max<std::size_t>(1024, n * 16ull);
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
        for (auto& t : temp) insertBody(root, &t, pool, poolMaxNodes);
        computeMass(root);

        double maxA = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            auto p = temp[i].GetPos();
            double ax = 0.0, ay = 0.0, az = 0.0;
            accumulateAccel(temp[i], root, p.x, p.y, p.z, ax, ay, az, theta, softSq);
            outAcc[i] = glm::dvec3(ax, ay, az);
            maxA = std::max(maxA, glm::length(outAcc[i]));
        }
        if (outMaxAcc) *outMaxAcc = maxA;
    };

    double remaining = dt;
    constexpr int MAX_SUBSTEPS = 1024;
    int substeps = 0;
    while (remaining > 0.0 && substeps < MAX_SUBSTEPS) {
        double maxAcc = 0.0;
        computeAccelerationsBH(x, k1v, &maxAcc);
        for (std::size_t i = 0; i < n; ++i) k1x[i] = v[i];

        double suggested = remaining;
        if (maxAcc > 0.0) {
            suggested = std::min(suggested, 0.02 / std::sqrt(maxAcc));
        }
        const double minStep = dt / 4096.0;
        double h = std::clamp(suggested, minStep, remaining);

        for (std::size_t i = 0; i < n; ++i) {
            tmpX[i] = x[i] + 0.5 * h * k1x[i];
            tmpV[i] = v[i] + 0.5 * h * k1v[i];
        }
        computeAccelerationsBH(tmpX, k2v, nullptr);
        for (std::size_t i = 0; i < n; ++i) k2x[i] = tmpV[i];

        for (std::size_t i = 0; i < n; ++i) {
            tmpX[i] = x[i] + 0.5 * h * k2x[i];
            tmpV[i] = v[i] + 0.5 * h * k2v[i];
        }
        computeAccelerationsBH(tmpX, k3v, nullptr);
        for (std::size_t i = 0; i < n; ++i) k3x[i] = tmpV[i];

        for (std::size_t i = 0; i < n; ++i) {
            tmpX[i] = x[i] + h * k3x[i];
            tmpV[i] = v[i] + h * k3v[i];
        }
        computeAccelerationsBH(tmpX, k4v, nullptr);
        for (std::size_t i = 0; i < n; ++i) k4x[i] = tmpV[i];

        const double w = h / 6.0;
        for (std::size_t i = 0; i < n; ++i) {
            x[i] += w * (k1x[i] + 2.0 * k2x[i] + 2.0 * k3x[i] + k4x[i]);
            v[i] += w * (k1v[i] + 2.0 * k2v[i] + 2.0 * k3v[i] + k4v[i]);
        }

        remaining -= h;
        ++substeps;
    }

    computeAccelerationsBH(x, k1v, nullptr);
    for (std::size_t i = 0; i < n; ++i) {
        objs[i].position = x[i];
        objs[i].velocity = v[i];
        objs[i].acceleration = k1v[i];
    }
}
