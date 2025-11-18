// barnes_hut_step.cpp
#include <vector>
#include <cmath>
#include <algorithm>
#include "object.hpp"
#include "barnes_hut.hpp"
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

static constexpr double G_KM = 6.6743e-20; // км^3 / (кг·с^2)
static constexpr double THETA = 0.5;

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

static void insertBody(OctreeNode* node, Object* body) {
    auto p = body->GetPos();
    if (!node->children[0] && !node->body) {
        node->body = body;
        return;
    }
    if (!node->children[0] && node->body) {
        Object* old = node->body;
        node->body = nullptr;
        double half = node->size * 0.5;
        for (int i = 0; i < 8; ++i) {
            double nx = node->x0 + ((i & 1) ? half : 0.0);
            double ny = node->y0 + ((i & 2) ? half : 0.0);
            double nz = node->z0 + ((i & 4) ? half : 0.0);
            node->children[i] = new OctreeNode(nx, ny, nz, half);
        }
        int o = getOctant(node, old->GetPos().x, old->GetPos().y, old->GetPos().z);
        node->children[o]->body = old;
    }
    int idx = getOctant(node, p.x, p.y, p.z);
    OctreeNode* child = node->children[idx];
    if (!child->children[0] && !child->body) {
        child->body = body;
    } else {
        insertBody(child, body);
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
                            double& ax, double& ay, double& az)
{
    if (node->mass <= 0.0) return;

    if (!node->children[0]) {
        if (node->body == &obj) return;
        double dx = node->cx - px;
        double dy = node->cy - py;
        double dz = node->cz - pz;
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 == 0.0) return;
        double invR = 1.0 / std::sqrt(r2);
        double invR3 = invR * invR * invR;
        double f = G_KM * node->mass * invR3;
        ax += f * dx;
        ay += f * dy;
        az += f * dz;
        return;
    }

    double dx = node->cx - px;
    double dy = node->cy - py;
    double dz = node->cz - pz;
    double r2 = dx*dx + dy*dy + dz*dz;
    if (r2 == 0.0) return;
    double r = std::sqrt(r2);

    if ((node->size / r) < THETA) {
        double invR = 1.0 / r;
        double invR3 = invR * invR * invR;
        double f = G_KM * node->mass * invR3;
        ax += f * dx;
        ay += f * dy;
        az += f * dz;
    } else {
        for (int i = 0; i < 8; ++i) {
            if (node->children[i]) {
                accumulateAccel(obj, node->children[i], px, py, pz, ax, ay, az);
            }
        }
    }
}

static void deleteTree(OctreeNode* node) {
    for (int i = 0; i < 8; ++i) {
        if (node->children[i]) deleteTree(node->children[i]);
    }
    delete node;
}

void simulationStepBarnesHutCPU(std::vector<Object>& objs, float dt, bool pause) {
    if (pause || objs.empty()) return;

    auto p0 = objs[0].GetPos();
    double minX = p0.x, maxX = p0.x;
    double minY = p0.y, maxY = p0.y;
    double minZ = p0.z, maxZ = p0.z;

    for (const auto& o : objs) {
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

    OctreeNode* root = new OctreeNode(cx - side * 0.5,
                                      cy - side * 0.5,
                                      cz - side * 0.5,
                                      side);

    for (auto& o : objs) insertBody(root, &o);
    computeMass(root);

    for (auto& o : objs) {
        auto p = o.GetPos();
        double ax = 0.0, ay = 0.0, az = 0.0;
        accumulateAccel(o, root, p.x, p.y, p.z, ax, ay, az);
        o.accelerate(static_cast<float>(ax),
                     static_cast<float>(ay),
                     static_cast<float>(az),
                     dt);
        o.UpdatePos(dt);
    }

    deleteTree(root);
}
