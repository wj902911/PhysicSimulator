#include "AABB.h"
#include "CollisionDetection.h"
#include "RigidBodyInstance.h"
#include "RigidBodyTemplate.h"
#include "../shared/VectorMath.h"
#include "../shared/rpoly.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream>
#include <unordered_set>

namespace cloth2 {

using namespace std;

ColInfo findCollisions(Eigen::Ref<const Eigen::MatrixX3d> oldV,
                       Eigen::Ref<const Eigen::MatrixX3d> V,
                       Eigen::Ref<const Eigen::MatrixX3i> F,
                       const std::vector<std::shared_ptr<RigidBodyInstance>>& rbodies)
{
    ColInfo ret;
    auto& vertCollisions = ret.vertCollisions;
    auto& faceCollisions = ret.faceCollisions;
    vertCollisions.resize(V.rows());
    faceCollisions.resize(F.rows());
    auto clothAABB = buildAABB(oldV, V, F);

    for (int body = 0; body < rbodies.size(); body++) {
        std::vector<Collision> potentialCollisions;
        // TODO: Implement getPBDInfo
        auto pbd_tup = rbodies[body]->getPBDInfo();
        const auto& oldPos = std::get<0>(pbd_tup);
        const auto& newPos = std::get<1>(pbd_tup);
        // TODO: Broad Phase, cloth vs body
        //       Construct Rigid body's AABB here
        auto bodyAABB = buildAABB(oldPos, newPos, rbodies[body]->getTemplate().getFaces());
        intersect(clothAABB, bodyAABB, potentialCollisions);
        const Eigen::MatrixX3i& bodyF = rbodies[body]->getTemplate().getFaces();
        for (const auto& it : potentialCollisions) {
            Eigen::Vector3i face = bodyF.row(it.collidingTriangle2);
            for (int vert = 0; vert < 3; vert++) {
                int vidx = F(it.collidingTriangle1, vert);
                double t;
                // TODO Implement Vertex-Face CTCD
                if (CTCD::vertexFaceCTCD(
                                         oldV.row(vidx).transpose(),
                                         oldPos.row(face(0)).transpose(),
                                         oldPos.row(face(1)).transpose(),
                                         oldPos.row(face(2)).transpose(),
                                         V.row(vidx).transpose(),
                                         newPos.row(face(0)).transpose(),
                                         newPos.row(face(1)).transpose(),
                                         newPos.row(face(2)).transpose(),
                                         1e-6,
                                         t)) {
                    vertCollisions[vidx].emplace(body, it.collidingTriangle2);
                    //std::cout << "collision detected!" << std::endl;
                }
            }
        }
        int facecnt = 0;
        // TODO: Broad Phase, body vs cloth
        //       Can reuse the constructed Rigid body's AABB.
        potentialCollisions.clear();
        intersect(bodyAABB, clothAABB, potentialCollisions);
        for (const auto& it : potentialCollisions) {
            Eigen::Vector3i face = F.row(it.collidingTriangle2);
            for (int vert = 0; vert < 3; vert++) {
                int vidx = bodyF(it.collidingTriangle1, vert);
                double t;
                // TODO Implement Vertex-Face CTCD
                if (CTCD::vertexFaceCTCD(
                                         oldPos.row(vidx).transpose(),
                                         oldV.row(face(0)).transpose(),
                                         oldV.row(face(1)).transpose(),
                                         oldV.row(face(2)).transpose(),
                                         newPos.row(vidx).transpose(),
                                         V.row(face(0)).transpose(),
                                         V.row(face(1)).transpose(),
                                         V.row(face(2)).transpose(),
                                         1e-6,
                                         t)) {
                    faceCollisions[it.collidingTriangle2].emplace(body, vidx);
                    //std::cout << "collision detected!" << std::endl;
                    facecnt++;
                }
            }
        }
    }
    return ret;
}


// TODO: Implmentation of Vertex-Face CTCD Narrow Phase
bool CTCD::vertexFaceCTCD(const Eigen::Vector3d& q0start,//a0
                          const Eigen::Vector3d& q1start,//b0
                          const Eigen::Vector3d& q2start,//c0
                          const Eigen::Vector3d& q3start,//d0
                          const Eigen::Vector3d& q0end,//a1
                          const Eigen::Vector3d& q1end,//b1
                          const Eigen::Vector3d& q2end,//c1
                          const Eigen::Vector3d& q3end,//d1
                          double eta,
                          double& t)
{
    Eigen::Vector4d coeffs;
    coeffs.setZero();
    coeffs += coeffCalc(q1start, q1end, q2start, q2end, q3start, q3end);
    coeffs -= coeffCalc(q1start, q1end, q0start, q0end, q3start, q3end);
    coeffs -= coeffCalc(q0start, q0end, q2start, q2end, q3start, q3end);
    coeffs -= coeffCalc(q1start, q1end, q2start, q2end, q0start, q0end);

    double coef[4] = { 0, 0, 0, 0 };
    int degree = 3;

    RootFinder RF;
    if (coeffs(0)!=0) {
        coef[0] = coeffs(0);
        coef[1] = coeffs(1);
        coef[2] = coeffs(2);
        coef[3] = coeffs(3);
    } else {
        if (coeffs(1) != 0) {
            coef[0] = coeffs(1);
            coef[1] = coeffs(2);
            coef[2] = coeffs(3);
            degree = 2;
        } else {
            coef[0] = coeffs(2);
            coef[1] = coeffs(3);
            degree = 1;
        }
    }
    double zeror[3] = { 0, 0, 0 }, zeroi[3] = { 0, 0, 0 };
    RF.rpoly(coef, 3, zeror, zeroi);
    std::unordered_set<double> ts;
    std::unordered_set<double> ts_true;
    for (int i = 0; i < 3; i++) {
        if (zeroi[i] == 0 && zeror[i] > 0 && zeror[i] < 1)
            ts.insert(zeror[i]);
    }
    if (ts.empty())
        return false;
    else
        for (const auto& it : ts) {
            Eigen::Vector3d p = (1.0 - it) * q0start + it * q0end;
            Eigen::Vector3d a = (1.0 - it) * q1start + it * q1end;
            Eigen::Vector3d b = (1.0 - it) * q2start + it * q2end;
            Eigen::Vector3d c = (1.0 - it) * q3start + it * q3end;
            if (!pointInTriangle(a, b, c, p))
                ts.erase(it);
        }
    if (ts.empty())
        return false;
    else
        return true;
}

Eigen::Vector4d CTCD::coeffCalc(Eigen::Vector3d x0, Eigen::Vector3d x1, Eigen::Vector3d y0, Eigen::Vector3d y1, Eigen::Vector3d z0, Eigen::Vector3d z1)
{
    Eigen::Vector4d result;
    double a0 = x0.cross(y0).dot(z0);
    double a1 = x0.cross(y0).dot(z1);
    double a2 = x0.cross(y1).dot(z0);
    double a3 = x0.cross(y1).dot(z1);
    double a4 = x1.cross(y0).dot(z0);
    double a5 = x1.cross(y0).dot(z1);
    double a6 = x1.cross(y1).dot(z0);
    double a7 = x1.cross(y1).dot(z1);
    result(0) = -a0 + a1 + a2 - a3 + a4 - a5 - a6 + a7;
    result(1) = 3 * a0 - 2 * a1 - 2 * a2 + a3 - 2 * a4 + a5 + a6;
    result(2) = -3 * a0 + a1 + a2 + a4;
    result(3) = a0;
    return result;
}

bool CTCD::sameSide(Eigen::Vector3d ls, Eigen::Vector3d le, Eigen::Vector3d p1, Eigen::Vector3d p2)
{
    Eigen::Vector3d c1 = (le - ls).cross(p1 - ls);
    Eigen::Vector3d c2 = (le - ls).cross(p2 - ls);
    if (c1.dot(c2) >= 0)
        return true;
    else
        return false;
}

bool CTCD::pointInTriangle(Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c, Eigen::Vector3d p)
{
    if (sameSide(a, b, c, p) && sameSide(b, c, a, p) && sameSide(c, a, b, p))
        return true;
    else
        return false;
}

}
