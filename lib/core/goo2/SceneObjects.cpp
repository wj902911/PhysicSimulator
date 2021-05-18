#include "SceneObjects.h"
#include <iostream>

namespace goo2 {

void Spring::processSpringForce(const Eigen::VectorXd& q,
                                Eigen::Ref<Eigen::VectorXd> F,
                                TripletArray& H) const
{
    using namespace Eigen;

    Vector2d x1 = q.segment<2>(2 * p1);
    Vector2d x2 = q.segment<2>(2 * p2);
    double dist = (x2 - x1).norm();
    Vector2d localF = stiffness * (dist - restlen) / dist * (x2 - x1);
    F.segment<2>(2 * p1) += localF;
    F.segment<2>(2 * p2) -= localF;

    Matrix2d I;
    I.setIdentity();
    Matrix2d localH = stiffness * (1.0 - restlen / dist) * I;
    localH += stiffness * restlen * (x2 - x1) * (x2 - x1).transpose() / dist / dist / dist;

    for (int j = 0; j < 2; j++) {
        for (int k = 0; k < 2; k++) {
            H.emplace_back(2 * p1 + j, 2 * p1 + k, localH.coeff(j, k));
            H.emplace_back(2 * p2 + j, 2 * p2 + k, localH.coeff(j, k));
            H.emplace_back(2 * p1 + j, 2 * p2 + k, -localH.coeff(j, k));
            H.emplace_back(2 * p2 + j, 2 * p1 + k, -localH.coeff(j, k));
        }
    }
}

void Spring::processDampingForce(const Eigen::VectorXd& q,
                                 const Eigen::VectorXd& qprev,
                                 const SimParameters& params,
                                 Eigen::Ref<Eigen::VectorXd> F,
                                 TripletArray& H) const
{
    using namespace Eigen;
    Vector2d x1 = q.segment<2>(2 * p1);
    Vector2d x2 = q.segment<2>(2 * p2);
    Vector2d x1prev = qprev.segment<2>(2 * p1);
    Vector2d x2prev = qprev.segment<2>(2 * p2);

    Vector2d relvel = (x2 - x2prev) / params.timeStep - (x1 - x1prev) / params.timeStep;
    Vector2d localF = params.dampingStiffness * relvel;
    F.segment<2>(2 * p1) += localF;
    F.segment<2>(2 * p2) -= localF;

    Matrix2d I;
    I.setIdentity();
    Matrix2d localH = params.dampingStiffness * I / params.timeStep;

    for (int j = 0; j < 2; j++) {
        for (int k = 0; k < 2; k++) {
            H.emplace_back(2 * p1 + j, 2 * p1 + k, localH.coeff(j, k));
            H.emplace_back(2 * p2 + j, 2 * p2 + k, localH.coeff(j, k));
            H.emplace_back(2 * p1 + j, 2 * p2 + k, -localH.coeff(j, k));
            H.emplace_back(2 * p2 + j, 2 * p1 + k, -localH.coeff(j, k));
        }
    }
}

void RigidRod::computePenaltyForce(const Eigen::VectorXd& q,
                                   double penaltyStiffness,
                                   Eigen::Ref<Eigen::VectorXd> F) const
{
    // TODO: Implmentation Penalty Force
    using namespace Eigen;

    Vector2d x1 = q.segment<2>(2 * p1);
    Vector2d x2 = q.segment<2>(2 * p2);
    double squaredDist = (x2 - x1).squaredNorm();
    Vector2d localF = 4 * penaltyStiffness * (squaredDist - length * length) * (x2 - x1);
    F.segment<2>(2 * p1) += localF;
    F.segment<2>(2 * p2) -= localF;
}

}
