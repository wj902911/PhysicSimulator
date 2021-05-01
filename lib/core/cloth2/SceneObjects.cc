#include "SceneObjects.h"
#include <iostream>

namespace cloth2 {

JetSource::JetSource(const Eigen::Vector3d& source, // source
                     const Eigen::Vector3d& dir, // Direction
                     double jetStrength,
                     double jetAtten,
                     uint32_t id)
    : src_(source),
      dir_(dir),
      k_(jetStrength),
      alpha_(jetAtten),
      id_(id)
{
}

void JetSource::blow(const Eigen::Ref<const Eigen::MatrixXd> V,
                     const Eigen::Ref<const Eigen::MatrixXi> F,
                     const Eigen::Ref<const Eigen::VectorXd> A,
                     const Eigen::Ref<const Eigen::MatrixXd> N,
                     const Eigen::Ref<const Eigen::MatrixXd> C,
                     Eigen::Ref<Eigen::MatrixXd> outForces) const
{
    if (!enabled_)
        return ;
    double coe = k_ / 3.0;
    for (int f = 0; f < F.rows(); f++) {
        Eigen::Vector3d force = coe * A(f) *
                                N.row(f).dot(dir_) *
                                std::pow((C.row(f).transpose() - src_).dot(dir_), alpha_) *
                                N.row(f);
        for (int i = 0; i < 3; i++)
            outForces.row(F(f, i)) += force;
    }
}

}
