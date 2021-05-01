#ifndef CLOTH2_SCENE_OBJECTS_H
#define CLOTH2_SCENE_OBJECTS_H

#include <Eigen/Core>
#include <stdint.h>

namespace cloth2 {

class JetSource {
public:
    JetSource(const Eigen::Vector3d& source, // source
              const Eigen::Vector3d& dir, // Direction
              double jetStrength,
              double jetAtten,
              uint32_t id);
    uint32_t id() const { return id_; }
    
    void updateConfiguration(const Eigen::Vector3d& source,
                             const Eigen::Vector3d& dir)
    {
        src_ = source;
        dir_ = dir;
    }
    void updateParameter(double jetStrength,
                         double jetAtten)
    {
        k_ = jetStrength;
        alpha_ = jetAtten;
    }
    
    void enable() { enabled_ = true; }
    void disable()  { enabled_ = false; }

    void blow(const Eigen::Ref<const Eigen::MatrixXd> V,
              const Eigen::Ref<const Eigen::MatrixXi> F,
              const Eigen::Ref<const Eigen::VectorXd> A,
              const Eigen::Ref<const Eigen::MatrixXd> N,
              const Eigen::Ref<const Eigen::MatrixXd> C,
              Eigen::Ref<Eigen::MatrixXd> outForces) const;
private:
    Eigen::Vector3d src_;
    Eigen::Vector3d dir_;
    double k_;
    double alpha_;
    uint32_t id_;
    bool enabled_ = true;
};

}

#endif
