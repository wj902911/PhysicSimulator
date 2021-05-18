#ifndef PSYM_CORE_GOO2_GOOCORE_H
#define PSYM_CORE_GOO2_GOOCORE_H

#include "../PhysicsCore.h"
#include "SceneObjects.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/StdVector>

#include <memory>
#include <stdint.h>
#include <tuple>
#include <unordered_map>

namespace goo2 {

class SimParameters;
using TripletArray = std::vector<Eigen::Triplet<double>>;

class GooCore : public PhysicsCore {
public:
    GooCore();
    ~GooCore();

    virtual void initSimulation() override;

    virtual bool simulateOneStep() override;

    std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd>
    getCurrentMesh() const override;

    std::shared_ptr<SimParameters> getPointerToSimParameters()
    {
        return params_;
    }

    /*
     * addParticle: add a particle at (x, y)
     *
     * Returns: a list of non-negative IDs which uniquely identify the newly
     *          added particles during the lifetime of the simulation (i.e.
     *          until initSimulation() is called)
     *
     * Note: A flex rod returns mutliple UIDs. First one is the non-inert
     *       particle, and the remaining are the inert particles.
     */
    Eigen::Matrix<uint32_t, Eigen::Dynamic, 1>
    addParticle(double x, double y);

    /*
     * Basic Spring UI
     * addSaw: add a saw at (x, y)
     */
    void addSaw(double x, double y);

    /*
     * queryParticle: locate the particle with specific UID
     *
     * Returns: particle.
     *          particle may be nullptr if the particle with given uid cannot
     *          be found.
     */
    std::shared_ptr<Particle>
    queryParticle(uint32_t uid);

    auto
    listAllParticles() const { return particles_; }

    /*
     * Mapping between UID <-> offset in configuration vector
     *
     * Note, do NOT multiple the offset by 2. For example, the returned
     * configuration index for the first and second particle should be 0 and
     * 1, not 0 and 2.
     *
     * UID and offset can be distinguished by their types.
     * UID is uint32_t, and offset is int
     */
    int mapParticleToConfigurationIndex(uint32_t) const;
    uint32_t mapConfigurationIndexToParticle(int) const;

    /*
     * queryConnector: locate spring object from given pair of UIDs
     * Returns: Connector.
     *          Connector may be nullptr if the Connector connecting particles with
     *          given UIDs cannot be found.
     * CAVEAT: This function was querySpring in M2, but now we are adopting a
     *         more general interface.
     */
    std::shared_ptr<const Connector>
    queryConnector(uint32_t uid1, uint32_t uid2) const;

    /****************************************************
     * Internal function called directly or indirectly by simulateOneStep, and
     * they will be tested by the grading script.
     ****************************************************/

    /*
     * Build configurations from particle properties
     * Returns: tuple (q, qprev, v)
     */
    std::tuple<Eigen::VectorXd, // current position (q)
               Eigen::VectorXd, // previous position (qprev)
               Eigen::VectorXd> // velocity (v)
    buildConfiguration() const;

    /*
     * Unbuild Configuration vectors to particle properties
     * Arguments:
     *      q: current position vector
     *      lambda: current lambda vector
     *      v: current velocity vector
     */
    void unbuildConfiguration(const Eigen::Ref<Eigen::VectorXd> q,
                              const Eigen::Ref<Eigen::VectorXd> lambda,
                              const Eigen::Ref<Eigen::VectorXd> v);

    /*
     * Mass Matrix (M) is a sparse matrix, but we usually use M^{-1} instead.
     * Returns
     *      H: Triplets of the inverse mass matrix.
     */
    TripletArray computeMassInverseCoeff() const;

    /*
     * Return a vector to store the force.
     * Sancheck: this vector should have the same dimensions and types as of
     * the current position vector (q).
     */
    Eigen::VectorXd createZeroForce() const;

    /*
     * Arguments:
     *      F: [out] gravity force vector.
     */
    void processGravityForce(Eigen::Ref<Eigen::VectorXd> F) const;

    /*
     * Compute the spring force and its Hessian
     *
     * Arguments:
     *      q: [in] configuration vector.
     *      F: [out] spring force vector.
     * Returns:
     *      H: Triplets of the Hessian (dF)
     */
    TripletArray processSpringForce(const Eigen::VectorXd& q,
                                    Eigen::Ref<Eigen::VectorXd> F) const;

    /*
     * Compute the viscous damping force and its Hessian
     *
     * Arguments:
     *      q: [in] configuration vector.
     *      F: [out] spring force vector.
     * Returns:
     *      H: Triplets of the Hessian (dF)
     */
    TripletArray processDampingForce(const Eigen::VectorXd& q,
                                     const Eigen::VectorXd& qprev,
                                     Eigen::Ref<Eigen::VectorXd> F) const;

    /*
     * Compute the floor force and its Hessian
     *
     * Arguments:
     *      q: [in] configuration vector.
     *      F: [out] spring force vector.
     * Returns:
     *      H: Triplets of the Hessian (dF)
     */
    TripletArray processFloorForce(const Eigen::VectorXd& q,
                                   const Eigen::VectorXd& qprev,
                                   Eigen::Ref<Eigen::VectorXd> F) const;

    /****************************************************
     * MILESTONE 2 NEW FUNCTIONS
     ****************************************************/

    /*
     *** Penalty Force Method
     */
    /*
     * Compute penalty force
     *
     * Arguments:
     *      q: [in] configuration vector.
     *      F: [out] penalty force vector.
     */
    void computePenaltyForce(const Eigen::VectorXd& q,
                             Eigen::Ref<Eigen::VectorXd> F) const;

    /*
     *** Step and Project
     */
    /*
     * Return a zero vector to store the value of Lagrangian function (f in
     * the Section 4.1.2)
     * Sancheck: this vector should have the same dimension as of the total
     *           number of constraints plus the degrees of freedom
     */
    Eigen::VectorXd createZeroLagrangian() const;

    std::vector<std::shared_ptr<Connector>>
    listAllRigidRods();

    /*
     * Compute f and df for the step and project method
     *
     * Arguments:
     *      q: [in] the initial guess of the constrained q
     *      qunconstrained: [in] the unconstrained q
     *      lambda: [in] the initial guess of lambda
     *      Minv: [in] the inverse mass matrix
     *      f: [out] f vector
     * Returns:
     *      df: Triplets of the sparse matrix [df].
     */
    TripletArray computeProjectionInfo(const Eigen::VectorXd& q,
                                       const Eigen::VectorXd& qunconstrained,
                                       const Eigen::VectorXd& lambda,
                                       const Eigen::SparseMatrix<double>& Minv,
                                       Eigen::Ref<Eigen::VectorXd> f) const;

    /*
     * Compute f and df for the step and project method
     *
     * Arguments:
     *      q: [in] the initial guess of the constrained q
     *      lambda: [in] the initial guess of lambda
     *      Minv: [in] the inverse mass matrix
     * Returns:
     *      newq: Projected q
     *      new lambda: the actual value of lambda
     */
    std::tuple<Eigen::VectorXd, // newq
               Eigen::VectorXd> // new lambda
    projectOntoConstraints(const Eigen::VectorXd& q,
                           const Eigen::VectorXd& lambda,
                           const Eigen::SparseMatrix<double>& Minv) const;

    /*
     *** Constrained Lagrangian
     */
    /*
     * Return a vector to store the value of constraint function (g)
     * Sancheck: this vector should have the same dimension as of the total
     *           number of rigid rods.
     */
    Eigen::VectorXd createZeroConstraints() const;

    /*
     * Compute g and dg for the constrained lagrangian method
     *
     * Arguments:
     *      q: [in] the initial guess of the constrained q
     *      g: [out] g vector
     * Returns:
     *      dg: Triplets of the sparse matrix [dg].
     */
    TripletArray computeConstraintAndGradient(const Eigen::VectorXd& q,
                                              Eigen::Ref<Eigen::VectorXd> g) const;


    /*
     * Use Newton's method to solve f(\lambda^{i+1}) = g(\ldots) = 0
     *
     * Arguments:
     *      q: [in] the value of q^i
     *      lambda: [in] the initial guess of lambda
     *      v: [in] the value of v^i
     * Returns:
     *      Df: [optional] The df in the LAST iteration to compute lambda^{i+1}.
     *          We are going to use this, in conjuction with
     *          params_->NewtonMaxIters, for partial credits of your
     *          Constrained Lagrangian method. You may leave if empty if
     *          you're confident about your result.
     */
    Eigen::SparseMatrix<double> // Df
    solveConstrainedLagrangian(Eigen::Ref<Eigen::VectorXd> q,
                               Eigen::Ref<Eigen::VectorXd> lambda,
                               Eigen::Ref<Eigen::VectorXd> v,
                               const Eigen::SparseMatrix<double>& Minv) const;

    /*
     *** Flexible Rods
     */

    auto
    listAllHinges() const { return bendingStencils_; }

    /*
     * Compute the effective mass of a particle
     *
     * Arguments:
     *      q: [in] Configuration index of the particle.
     * Returns:
     *      The corresponding particle's effective mass
     */
    double getTotalParticleMass(int idx) const;

    /*
     * Compute the bending force
     *
     * Arguments:
     *      q: [in] configuration vector.
     *      F: [out] bending force vector.
     */
    void processBendingForce(const Eigen::VectorXd& q,
                             Eigen::Ref<Eigen::VectorXd> F) const;

private:
    /****************************************************
     * Internal function called directly or indirectly by simulateOneStep, but
     * we are not going to test them.
     ****************************************************/
    void computeMassInverse(Eigen::SparseMatrix<double>& Minv);
    void numericalIntegration(Eigen::VectorXd& q,
                              Eigen::VectorXd& lambda, // Note: this takes lambda instead of qprev
                              Eigen::VectorXd& v);

    void computeForceWithoutHessian(const Eigen::VectorXd& q,
                                    const Eigen::VectorXd& qprev,
                                    Eigen::VectorXd& F) const;

    double ptSegmentDist(const Eigen::Vector2d& p,
                         const Eigen::Vector2d& q1,
                         const Eigen::Vector2d& q2) const;

    void detectSawedConnectors(std::set<int>& connectorsToDelete);
    void detectSawedParticles(std::set<int>& particlesToDelete);
    void deleteSawedObjects();
    void pruneOverstrainedSprings();

private:
    uint32_t particle_unique_id_;
    std::shared_ptr<SimParameters> params_;
    double time_;
    // This is required if we want mutable particles in python
    std::vector<std::shared_ptr<Particle>> particles_;
    std::vector<std::shared_ptr<Connector>> connectors_;
    std::vector<Saw> saws_;
    std::vector<BendingStencil> bendingStencils_;

    std::unordered_map<uint32_t, int> particle_offset_;

    int getNumRigidRods() const;
};

}

#endif
