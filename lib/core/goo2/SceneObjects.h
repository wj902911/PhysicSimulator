#ifndef PSYM_CORE_GOO2_SCENEOBJECTS_H
#define PSYM_CORE_GOO2_SCENEOBJECTS_H

#include "SimParameters.h"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <limits>
#include <cmath>
#include <set>
#include <vector>

namespace goo2 {

using TripletArray = std::vector<Eigen::Triplet<double>>;

struct Particle {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Particle(Eigen::Vector2d pos, double mass, bool isFixed, bool isInert, uint32_t uid)
        : pos(pos)
        , mass(mass)
        , fixed(isFixed)
        , inert(isInert)
        , uid(uid)
    {
        vel.setZero();
    }

    Eigen::Vector2d pos;
    Eigen::Vector2d vel;
    double mass;
    bool fixed;
    bool inert;

    uint32_t uid = -1;
};

struct Connector {
public:
    Connector(int p1, int p2, double mass)
        : p1(p1)
        , p2(p2)
        , mass(mass)
    {
    }
    virtual ~Connector() {}

    virtual SimParameters::ConnectorType getType() const = 0;

    int p1;
    int p2;
    double mass;

    std::set<int> associatedBendingStencils;

    /*
     * Per-connector spring force.
     *
     * Arguments:
     *      q: configuration vector of current positions
     *      F: [out] force vector
     *      H: [out] Hessian triplets
     *
     * You are supposed to override this function at Spring class
     */
    virtual void processSpringForce(const Eigen::VectorXd& q,
                                    Eigen::Ref<Eigen::VectorXd> F,
                                    TripletArray& H) const {}
    /*
     * Per-connector dampling force
     *
     * Arguments:
     *      q: configuration vector of current positions
     *      qprev: configuration vector of previous positions
     *      params: simulation parameter
     *      F: [out] force vector
     *      H: [out] Hessian triplets
     *
     * You are supposed to override this function at Spring class
     */
    virtual void processDampingForce(const Eigen::VectorXd& q,
                                     const Eigen::VectorXd& qprev,
                                     const SimParameters& params,
                                     Eigen::Ref<Eigen::VectorXd> F,
                                     TripletArray& H) const {}

    /*
     * Per-connector penalty force.
     *
     * Arguments:
     *      q: configuration vector of current positions
     *      F: [out] force vector
     *
     * You are supposed to override this function at RigidRod class
     */
    virtual void computePenaltyForce(const Eigen::VectorXd& q,
                                      double penaltyStiffness,
                                      Eigen::Ref<Eigen::VectorXd> F) const
    {
    }

    /*
     * Load the lambda
     *
     * NAN will be returned if the corresponding connector doesn't have lambda
     * associated with.
     */
    virtual double loadLambda() const { return std::numeric_limits<double>::quiet_NaN(); }

    /*
     * Store the lambda to the connector
     *
     * Arguments:
     *      lambdas: [in] the vector of lambdas
     *      i: [in] the index of the value in lambdas vector.
     * Returns:
     *      The number of lambdas stored in the connector.
     *      0 shall be returned if the corresponding connector doesn't have
     *      lambda associated with.
     *      TIPS: you only need to return 1 in some subclass of Connector in this project.
     */
    virtual int storeLambda(const Eigen::Ref<const Eigen::VectorXd> lambdas, int ) { return 0; }
};

struct Spring : public Connector {
public:
    Spring(int p1, int p2, double mass, double stiffness, double restlen, bool canSnap)
        : Connector(p1, p2, mass)
        , stiffness(stiffness)
        , restlen(restlen)
        , canSnap(canSnap)
    {
    }

    virtual SimParameters::ConnectorType getType() const override
    {
        return SimParameters::CT_SPRING;
    }

    double stiffness;
    double restlen;
    bool canSnap;

    void processSpringForce(const Eigen::VectorXd& q,
                            Eigen::Ref<Eigen::VectorXd> F,
                            TripletArray& H) const override;
    void processDampingForce(const Eigen::VectorXd& q,
                             const Eigen::VectorXd& qprev,
                             const SimParameters& params,
                             Eigen::Ref<Eigen::VectorXd> F,
                             TripletArray& H) const override;
};

struct RigidRod : public Connector {
public:
    RigidRod(int p1, int p2, double mass, double length)
        : Connector(p1, p2, mass)
        , lambda(0)
        , length(length)
    {
    }

    virtual SimParameters::ConnectorType getType() const override
    {
        return SimParameters::CT_RIGIDROD;
    }

    virtual void computePenaltyForce(const Eigen::VectorXd& q,
                                     double penaltyStiffness,
                                     Eigen::Ref<Eigen::VectorXd> F) const override;

    virtual double loadLambda() const override { return lambda; }
    virtual int storeLambda(const Eigen::Ref<const Eigen::VectorXd> lambdas, int index) override
    {
        this->lambda = lambdas(index);
        return 1;
    }

    double lambda;
    double length;
};

struct Saw {
public:
    Saw(Eigen::Vector2d pos, double radius)
        : pos(pos)
        , radius(radius)
    {
    }

    Eigen::Vector2d pos;
    double radius;
};

struct BendingStencil {
public:
    BendingStencil(int p1, int p2, int p3, double kb)
        : p1(p1)
        , p2(p2)
        , p3(p3)
        , kb(kb)
    {
    }

    int p1, p2, p3;
    double kb;
};

}

#endif
