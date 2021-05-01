#ifndef PSIM_CORE_CLOTH2_RIGIDBODYINSTANCE_H
#define PSIM_CORE_CLOTH2_RIGIDBODYINSTANCE_H

#include <Eigen/Core>
#include <list>
#include <vector>
#include <tuple>
#include <memory>

namespace cloth2 {

class RigidBodyTemplate;
class AABBNode;
struct ColInfo;

class RigidBodyInstance
{
public:
    RigidBodyInstance(std::shared_ptr<RigidBodyTemplate> rbtemplate,
                      const Eigen::Vector3d &c,
		      const Eigen::Vector3d &theta,
		      const Eigen::Vector3d &cvel,
		      const Eigen::Vector3d &w,
		      double density);
    ~RigidBodyInstance();

    Eigen::Vector3d c;
    Eigen::Vector3d theta;

    Eigen::Vector3d cvel;
    Eigen::Vector3d w;

    double density;
    
    const RigidBodyTemplate &getTemplate() const {return *rbtemplate_;}
    std::shared_ptr<RigidBodyTemplate> getTemplatePointer() const {return rbtemplate_;}
    
    Eigen::Vector3d oldc;
    Eigen::Vector3d oldtheta;
    
    std::tuple<Eigen::MatrixX3d, // old position
               Eigen::MatrixX3d> // new position
    getPBDInfo() const;
    std::tuple<Eigen::Vector3d,
               Eigen::Vector3d>
    getVertPos(int vertIndex) const;

private:
    std::shared_ptr<RigidBodyTemplate> rbtemplate_;
};

}

#endif // RIGIDBODYINSTANCE_H
