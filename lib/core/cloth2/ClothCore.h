#ifndef PSIM_CORE_CLOTH2_CLOTHCORE_H
#define PSIM_CORE_CLOTH2_CLOTHCORE_H

#include "../PhysicsCore.h"
#include <Eigen/SparseCore>
#include <memory>
#include <vector>
#include <tuple>
#include <stdint.h>
#include "RigidBodyInstance.h"
#include "RigidBodyTemplate.h"

namespace cloth2 {

struct SimParameters;
class JetSource;
struct ColInfo;

struct IntermediateRecord;

class ClothCore : public PhysicsCore {
public:
    ClothCore();
    ~ClothCore();

    void attachMesh(Eigen::Ref<Eigen::MatrixX3d> V,
                    Eigen::Ref<Eigen::MatrixX3i> F,
                    double scale = 50.0);

    void attachBody(Eigen::Ref<Eigen::MatrixX3d> V,
                    Eigen::Ref<Eigen::MatrixX3i> F,
                    double scale = 1.0);

    virtual void initSimulation() override;

    std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd>
    getCurrentMesh() const override;

    void
    setSystemPhase(Eigen::Ref<const Eigen::MatrixXd> Q,
                   Eigen::Ref<const Eigen::MatrixXd> Qdot)
    {
        Q_ = Q;
        Qdot_ = Qdot;
    }

    std::tuple<Eigen::MatrixXd, // vertices
               Eigen::MatrixXd> // velocities
    getSystemPhase() const
    {
        return std::make_tuple(Q_, Qdot_);
    }

    Eigen::MatrixXi getHinges() const
    {
        return H_;
    }

    void computeHinges();
    void applyConstraints();
    void applyPinConstraints();
    void applyStretchConstraints();
    void applyBendingConstraints();
    void applyDragConstraints(int vertex);

    uint32_t addJetSource(const Eigen::Vector3d& xc, // Source
                          const Eigen::Vector3d& dir, // Direction
                          double jetStrength,
                          double jetAtten);
    std::shared_ptr<JetSource> queryJetSource(uint32_t);
    void removeJetSource(uint32_t);
    void clearJetSource();

    Eigen::MatrixXd computeJetForces() const;

    std::tuple<Eigen::Vector3d, // oldcentroid
               Eigen::Vector3d, // newcentroid
               Eigen::Matrix3d> // The solution
    solveOrthogonalProcrustes(Eigen::VectorXi face_id) const;

    std::shared_ptr<SimParameters> getPointerToSimParameters()
    {
        return params_;
    }

    virtual bool simulateOneStep() override;

    void hold(int vertex_id, const Eigen::Vector3d& position);
    void updateHold(const Eigen::Vector3d& position);
    void releaseHold();

    void cutFaces(const std::vector<int>& faces);

    void deployBevo();
    void launchBevo(Eigen::Vector3d src,
                    Eigen::Vector3d dir,
                    Eigen::Vector3d orientation,
                    Eigen::Vector3d omega);

    void applyClothVertexVsBodyFaceConstraints(const ColInfo& col_info);
    void applyClothFaceVsBodyVertexConstraints(const ColInfo& col_info);

    std::vector<std::shared_ptr<RigidBodyInstance>> listRigidBodyInstances() const {
        return bodyInstances_;
    }

    Eigen::MatrixXd record_oldQ;
    std::vector<IntermediateRecord> records;
private:
    std::shared_ptr<SimParameters> params_;

    Eigen::MatrixXd origQ_;
    Eigen::MatrixXd Q_;
    Eigen::MatrixXd Qdot_;
    Eigen::MatrixXi F_;
    Eigen::MatrixXi origF_;

    Eigen::MatrixXi H_;

    std::vector<int> pinnedVerts_;

    int clickedVertex_;
    Eigen::Vector3d curPos_;
    uint32_t jets_id_ = 0;
    std::vector<std::shared_ptr<JetSource>> jets_;

    std::shared_ptr<RigidBodyTemplate> body_;
    std::shared_ptr<RigidBodyTemplate> launchBody_;
    std::vector<std::shared_ptr<RigidBodyInstance>> bodyInstances_;

    Eigen::SparseMatrix<double> Minv_;

    bool bevoDeployed_;
};

}

#endif
