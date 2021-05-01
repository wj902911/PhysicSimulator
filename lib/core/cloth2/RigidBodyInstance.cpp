#include "RigidBodyInstance.h"
#include "../shared/VectorMath.h"
#include "RigidBodyTemplate.h"
#include "AABB.h"
#include "CollisionDetection.h"
#include <Eigen/Geometry>
#include <iostream>
//#include "core/shared/VectorMath.h"

using namespace Eigen;
using namespace std;

namespace cloth2 {

RigidBodyInstance::RigidBodyInstance(std::shared_ptr<RigidBodyTemplate> rbtemplate,
    const Eigen::Vector3d &c, const Eigen::Vector3d &theta,
    const Eigen::Vector3d &cvel, const Eigen::Vector3d &w,
    double density)
    : c(c), theta(theta), cvel(cvel), w(w), density(density), rbtemplate_(rbtemplate)
{
}

RigidBodyInstance::~RigidBodyInstance()
{    
}

std::tuple<Eigen::MatrixX3d, // old position
           Eigen::MatrixX3d> // new position
RigidBodyInstance::getPBDInfo() const
{
    int nverts = getTemplate().getVerts().rows();
    Eigen::MatrixXd oldPos(nverts, 3);
    Eigen::MatrixXd newPos(nverts, 3);
    // TODO: Implmentation of additional member functions
#if 1
    Eigen::Matrix3d I;
    I.setIdentity();
    Eigen::Matrix3d oldThetaX = VectorMath::crossProductMatrix(oldtheta.normalized());
    Eigen::Matrix3d thetaX = VectorMath::crossProductMatrix(theta.normalized());
    double cosOldThetaN = cos(oldtheta.norm());
    double cosThetaN = cos(theta.norm());
    double sinOldThetaN = sin(oldtheta.norm());
    double sinThetaN = sin(theta.norm());
    Eigen::Matrix3d rotOldTheta = cosOldThetaN * I + sinOldThetaN * oldThetaX + (1 - cosOldThetaN) * oldtheta.normalized() * oldtheta.normalized().transpose();
    Eigen::Matrix3d rotTheta = cosThetaN * I + sinThetaN * thetaX + (1 - cosThetaN) * theta.normalized() * theta.normalized().transpose();
    for (int i = 0; i < nverts; i++) {
        oldPos.row(i) = (oldc + rotOldTheta * getTemplate().getVerts().row(i).transpose()).transpose();
        newPos.row(i) = (c + rotTheta * getTemplate().getVerts().row(i).transpose()).transpose();
    }
#endif
#if 0
    for (int i = 0; i < nverts; i++)
        std::tie(oldPos.row(i), newPos.row(i)) = getVertPos(i);
#endif
    return std::make_tuple(oldPos, newPos);
}

std::tuple<Eigen::Vector3d, Eigen::Vector3d> RigidBodyInstance::getVertPos(int vertIndex) const
{
    Eigen::Vector3d oldPos;
    Eigen::Vector3d newPos;
    Eigen::Matrix3d I;
    I.setIdentity();
    Eigen::Matrix3d oldThetaX = VectorMath::crossProductMatrix(oldtheta.normalized());
    Eigen::Matrix3d thetaX = VectorMath::crossProductMatrix(theta.normalized());
    double cosOldThetaN = cos(oldtheta.norm());
    double cosThetaN = cos(theta.norm());
    double sinOldThetaN = sin(oldtheta.norm());
    double sinThetaN = sin(theta.norm());
    Eigen::Matrix3d rotOldTheta = cosOldThetaN * I + sinOldThetaN * oldThetaX + (1 - cosOldThetaN) * oldtheta.normalized() * oldtheta.normalized().transpose();
    Eigen::Matrix3d rotTheta = cosThetaN * I + sinThetaN * thetaX + (1 - cosThetaN) * theta.normalized() * theta.normalized().transpose();
    oldPos = (oldc + rotOldTheta * getTemplate().getVerts().row(vertIndex).transpose()).transpose();
    newPos = (c + rotTheta * getTemplate().getVerts().row(vertIndex).transpose()).transpose();

    return std::make_tuple(oldPos, newPos);
}

}
