#ifndef PSIM_CORE_CLOTH2_COLLISION_DETECTION_H
#define PSIM_CORE_CLOTH2_COLLISION_DETECTION_H

#include <Eigen/Core>
#include <set>
#include <vector>
#include <memory>

namespace cloth2 {

class RigidBodyInstance;
struct AABBNode;

struct ColInfo {
    std::vector<std::set<std::pair<int, int>>> vertCollisions;
    std::vector<std::set<std::pair<int, int>>> faceCollisions;

    void resize(size_t s) {
        vertCollisions.resize(s);
        faceCollisions.resize(s);
    }

    void vert_add(size_t off, int first, int second) {
        vertCollisions[off].emplace(first, second);
    }

    void face_add(size_t off, int first, int second) {
        faceCollisions[off].emplace(first, second);
    }
};

ColInfo findCollisions(Eigen::Ref<const Eigen::MatrixX3d> oldV,
                       Eigen::Ref<const Eigen::MatrixX3d> V,
                       Eigen::Ref<const Eigen::MatrixX3i> F,
                       const std::vector<std::shared_ptr<RigidBodyInstance>>& rbodies);

class CTCD {
public:
    // Vertex Face continunous collision detection
    // q0: vertex
    // q1-q3: faces
    //
    // For partial credit, the following inputs are tested
    //  1. q0start == q0end, or
    //  2. q1/2/3start == q1/2/3end
    // In other words, either the vertex or the face is stationary.
    static bool vertexFaceCTCD(const Eigen::Vector3d& q0start,
                               const Eigen::Vector3d& q1start,
                               const Eigen::Vector3d& q2start,
                               const Eigen::Vector3d& q3start,
                               const Eigen::Vector3d& q0end,
                               const Eigen::Vector3d& q1end,
                               const Eigen::Vector3d& q2end,
                               const Eigen::Vector3d& q3end,
                               double eta,
                               double& t);

    static Eigen::Vector4d coeffCalc(Eigen::Vector3d x0, Eigen::Vector3d x1, Eigen::Vector3d y0, Eigen::Vector3d y1, Eigen::Vector3d z0, Eigen::Vector3d z1);

    static bool sameSide(Eigen::Vector3d ls, Eigen::Vector3d le, Eigen::Vector3d p1, Eigen::Vector3d p2);

    static bool pointInTriangle(Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c, Eigen::Vector3d p);
};

}

#endif // PSIM_CORE_CLOTH2_COLLISION_DETECTION_H
