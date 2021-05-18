#include "GooCore.h"
#include "SimParameters.h"
#include <Eigen/Geometry>
#include <Eigen/SparseQR>
#include <algorithm>
#include <unordered_map>
#include <iostream>

using namespace Eigen;

namespace goo2 {

SparseMatrix<double>& operator<<(SparseMatrix<double>& sp, const TripletArray& tris)
{
    sp.setFromTriplets(tris.begin(), tris.end());
    return sp;
}

GooCore::GooCore()
{
    params_ = std::make_shared<SimParameters>();
}

GooCore::~GooCore()
{
}

void GooCore::initSimulation()
{
    particle_unique_id_ = 0;
    time_ = 0;
    particles_.clear();
    connectors_.clear();
    saws_.clear();
    bendingStencils_.clear();
    // TODO: other initializatios
}

bool GooCore::simulateOneStep()
{
    VectorXd q, lambda, v;
    std::tie(q, lambda, v) = buildConfiguration();
    numericalIntegration(q, lambda, v);
    unbuildConfiguration(q, lambda, v);

    pruneOverstrainedSprings();
    deleteSawedObjects();
    time_ += params_->timeStep;
    return false;
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd>
GooCore::getCurrentMesh() const
{
    double baseradius = 0.02;
    double pulsefactor = 0.1;
    double pulsespeed = 50.0;

    int sawteeth = 20;
    double sawdepth = 0.1;
    double sawangspeed = 10.0;

    double baselinewidth = 0.005;

    int numcirclewedges = 20;

    // this is terrible. But, easiest to get up and running

    std::vector<Eigen::Vector3d> verts;
    std::vector<Eigen::Vector3d> vertexColors;
    std::vector<Eigen::Vector3i> faces;

    int idx = 0;

    double eps = 1e-4;

    if (params_->floorEnabled) {
        for (int i = 0; i < 6; i++) {
            vertexColors.emplace_back(0.3, 1.0, 0.3);
        }

        verts.emplace_back(-1, -0.5, eps);
        verts.emplace_back(1, -0.5, eps);
        verts.emplace_back(-1, -1, eps);

        faces.emplace_back(idx, idx + 1, idx + 2);

        verts.emplace_back(-1, -1, eps);
        verts.emplace_back(1, -0.5, eps);
        verts.emplace_back(1, -1, eps);
        faces.emplace_back(idx + 3, idx + 4, idx + 5);
        idx += 6;
    }

    for (const auto& c : connectors_) {
        Vector2d sourcepos = particles_[c->p1]->pos;
        Vector2d destpos = particles_[c->p2]->pos;

        Vector2d vec = destpos - sourcepos;
        Vector2d perp(-vec[1], vec[0]);
        perp /= perp.norm();

        double dist = (sourcepos - destpos).norm();

        Eigen::Vector3d color;
        double width;
        switch (c->getType()) {
        case SimParameters::CT_SPRING: {
            if (c->associatedBendingStencils.empty())
                color << 0.0, 0.0, 1.0;
            else
                color << 0.75, 0.5, 0.75;
            width = baselinewidth / (1.0 + 20.0 * dist * dist);

            break;
        }
        case SimParameters::CT_RIGIDROD: {
            if (c->associatedBendingStencils.empty())
                color << 1.0, 0.0, 1.0;
            else
                color << 1.0, 1.0, 0.3;
            width = baselinewidth;
            break;
        }
        default:
            break;
        }

        for (int i = 0; i < 4; i++)
            vertexColors.emplace_back(color);

        verts.emplace_back(sourcepos[0] + width * perp[0], sourcepos[1] + width * perp[1], -eps);
        verts.emplace_back(sourcepos[0] - width * perp[0], sourcepos[1] - width * perp[1], -eps);
        verts.emplace_back(destpos[0] + width * perp[0], destpos[1] + width * perp[1], -eps);
        verts.emplace_back(destpos[0] - width * perp[0], destpos[1] - width * perp[1], -eps);

        faces.emplace_back(idx, idx + 1, idx + 2);
        faces.emplace_back(idx + 2, idx + 1, idx + 3);
        idx += 4;
    }

    int nparticles = particles_.size();

    for (int i = 0; i < nparticles; i++) {
        double radius = baseradius * sqrt(getTotalParticleMass(i));
        radius *= (1.0 + pulsefactor * sin(pulsespeed * time_));

        Eigen::Vector3d color(0, 0, 0);

        if (particles_[i]->fixed) {
            radius = baseradius;
            color << 1.0, 0, 0;
        }

        for (int j = 0; j < numcirclewedges + 2; j++) {
            vertexColors.emplace_back(color);
        }

        verts.emplace_back(particles_[i]->pos[0], particles_[i]->pos[1], 0);

        const double PI = 3.1415926535898;
        for (int j = 0; j <= numcirclewedges; j++) {
            verts.emplace_back(particles_[i]->pos[0] + radius * cos(2 * PI * j / numcirclewedges),
                               particles_[i]->pos[1] + radius * sin(2 * PI * j / numcirclewedges),
                               0);
        }

        for (int j = 0; j <= numcirclewedges; j++) {
            faces.emplace_back(idx, idx + j + 1, idx + 1 + ((j + 1) % (numcirclewedges + 1)));
        }

        idx += numcirclewedges + 2;
    }

    for (const auto& saw : saws_) {
        double outerradius = saw.radius;
        double innerradius = (1.0 - sawdepth) * outerradius;

        Eigen::Vector3d color(0.5, 0.5, 0.5);

        int spokes = 2 * sawteeth;
        for (int j = 0; j < spokes + 2; j++) {
            vertexColors.emplace_back(color);
        }

        verts.emplace_back(saw.pos[0], saw.pos[1], 0);

        const double PI = 3.1415926535898;
        for (int i = 0; i <= spokes; i++) {
            double radius = (i % 2 == 0) ? innerradius : outerradius;
            verts.emplace_back(saw.pos[0] + radius * cos(2 * PI * i / spokes + sawangspeed * time_),
                               saw.pos[1] + radius * sin(2 * PI * i / spokes + sawangspeed * time_),
                               0.0);
        }

        for (int j = 0; j <= spokes; j++) {
            faces.emplace_back(idx, idx + j + 1, idx + 1 + ((j + 1) % (spokes + 1)));
        }

        idx += spokes + 2;
    }

    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;
    Eigen::MatrixXd renderC;

    renderQ.resize(verts.size(), 3);
    renderC.resize(vertexColors.size(), 3);
    for (int i = 0; i < verts.size(); i++) {
        renderQ.row(i) = verts[i];
        renderC.row(i) = vertexColors[i];
    }
    renderF.resize(faces.size(), 3);
    for (int i = 0; i < faces.size(); i++)
        renderF.row(i) = faces[i];
    return std::make_tuple(renderQ, renderF, renderC);
}

Eigen::Matrix<uint32_t, Eigen::Dynamic, 1>
GooCore::addParticle(double x, double y)
{
    Eigen::Matrix<uint32_t, Eigen::Dynamic, 1> ret;
    std::vector<int> ids;
    Vector2d newpos(x, y);
    double mass = params_->particleMass;
    if (params_->particleFixed)
        mass = std::numeric_limits<double>::infinity();

    int newid = particles_.size();
    particle_offset_[particle_unique_id_] = (int)particles_.size();
    particles_.emplace_back(std::make_shared<Particle>(newpos, mass, params_->particleFixed, false, particle_unique_id_));
    ids.emplace_back(particle_unique_id_);
    particle_unique_id_++;

    int numparticles = particles_.size() - 1;

    for (int i = 0; i < numparticles; i++) {
        if (particles_[i]->inert)
            continue;
        Vector2d pos = particles_[i]->pos;
        double dist = (pos - newpos).norm();
        if (dist > params_->maxSpringDist)
            continue;
        switch (params_->connectorType) {
            case SimParameters::CT_SPRING: {
                connectors_.emplace_back(std::make_shared<Spring>(newid, i, 0, params_->springStiffness / dist, dist, true));
                break;
            }
    // TODO: Handle Rigid rods and flexible rods.
            case SimParameters::CT_RIGIDROD: {
                connectors_.emplace_back(std::make_shared<RigidRod>(newid, i, 0, dist));
                break;
            }
            case SimParameters::CT_FLEXROD: {
                int nSecs = std::max(2, params_->rodSegments);
                double dist_j = dist / nSecs;
                double kb = params_->rodBendingStiffness / dist_j;
                Vector2d deltaPos = (pos - newpos).normalized() * dist_j;
                Vector2d pos_j = newpos;
                VectorXi pj(nSecs + 1);
                pj[0] = newid;
                pj[1] = particles_.size();
                for (int j = 2; j < nSecs; j++) {
                    pj[j] = pj[j - 1] + 1;
                }
                pj[nSecs] = i;
                //std::cout << pj.transpose() << std::endl;
                for (int j = 0; j < nSecs; j++) {
                    connectors_.emplace_back(std::make_shared<Spring>(pj[j], pj[j + 1], params_->rodDensity * dist_j, params_->rodStretchingStiffness / dist_j, dist_j, false));
                    //std::cout << pj[j] << "," << pj[j + 1] << std::endl;
                }
                int nConnecs = connectors_.size();
                for (int j = 0; j < nSecs - 1; j++) {
                    pos_j += deltaPos;
                    particle_offset_[particle_unique_id_] = (int)particles_.size();
                    particles_.emplace_back(std::make_shared<Particle>(pos_j, 0, false, true, particle_unique_id_));
                    //std::cout <<"inert particle added  "<< pos_j.transpose() << std::endl;
                    ids.emplace_back(particle_unique_id_);
                    particle_unique_id_++;
                    int nBendStens = bendingStencils_.size();
                    bendingStencils_.emplace_back(BendingStencil(pj[j], pj[j + 1], pj[j + 2], kb));
                    connectors_[nConnecs - nSecs + j]->associatedBendingStencils.insert(nBendStens);
                    connectors_[nConnecs - nSecs + j + 1]->associatedBendingStencils.insert(nBendStens);
                }
                break;
            }
            default:
                break;
        }
    }
    /*for (int i = 0; i < particles_.size(); i++) {
        std::cout << "pos:" << particles_[i]->pos.transpose() << ","
                  << "mass:" << particles_[i]->mass << ","
                  << "fixed:" << particles_[i]->fixed << ","
                  << "inert:" << particles_[i]->inert << ","
                  << "uid:" << particles_[i]->uid << std::endl;
    }*/
    for (const auto& pc : connectors_) {
        switch (pc->getType()) {
        case SimParameters::CT_SPRING: {
            std::shared_ptr<const Spring> ps = std::dynamic_pointer_cast<const Spring>(pc);
            std::cout << "Type:Spring, "
                      << "pos:" << ps->p1 << " " << ps->p2 << ", "
                      << "mass:" << ps->mass << ", "
                      << "stiffness:" << ps->stiffness << ", "
                      << "restlen:" << ps->restlen << ", "
                      << "canSnap:" << ps->canSnap << ", "
                      << "associatedBendingStencils:";
            for (std::set<int>::iterator it = pc->associatedBendingStencils.begin(); it != pc->associatedBendingStencils.end(); ++it) {
                std::cout << " " << *it;
            }
            std::cout << std::endl;
            break;
        }
        case SimParameters::CT_RIGIDROD: {
            std::shared_ptr<const RigidRod> pr = std::dynamic_pointer_cast<const RigidRod>(pc);
            std::cout << "Type:RigidRod, "
                      << "pos:" << pr->p1 << " " << pr->p2 << ", "
                      << "mass:" << pr->mass << ", "
                      << "restlen:" << pr->length << std::endl;
            break;
        }
        default:
            break;
        }
    }
    for (const auto& pb : bendingStencils_) {
        std::cout << "p1:" << pb.p1 << ","
                  << "p2:" << pb.p2 << ","
                  << "p3:" << pb.p3 << ","
                  << "kb:" << pb.kb << std::endl;
    }
    std::cout << std::endl;
    ret.resize(ids.size());
    for (int i = 0; i < ret.size(); i++)
        ret(i) = ids[i];
    return ret;
}

void GooCore::addSaw(double x, double y)
{
    saws_.emplace_back(Vector2d(x, y), params_->sawRadius);
}

std::shared_ptr<Particle>
GooCore::queryParticle(uint32_t uid)
{
    for (const auto& p : particles_)
        if (p->uid == uid)
            return p;
    return std::shared_ptr<Particle>(nullptr);
}

int GooCore::mapParticleToConfigurationIndex(uint32_t uid) const
{
    // TODO
    auto iter = particle_offset_.find(uid);
    //int p = iter->second;
    return (iter == particle_offset_.end()) ? -1 : iter->second;
}

uint32_t GooCore::mapConfigurationIndexToParticle(int i) const
{
    // TODO
    for (auto iter = particle_offset_.begin(); iter != particle_offset_.end(); iter++)
        if (iter->second == i)
            return iter->first;
    return 0;
}

std::shared_ptr<const Connector>
GooCore::queryConnector(uint32_t uid1, uint32_t uid2) const
{
    auto p1 = mapParticleToConfigurationIndex(uid1);
    auto p2 = mapParticleToConfigurationIndex(uid2);
    for (const auto& pc : connectors_) {
        const Connector& c = *pc;
        if (c.p1 == p1 && c.p2 == p2) {
            return pc;
        } else if (c.p2 == p1 && c.p1 == p2) {
            return pc;
        }
    }
    return std::shared_ptr<const Connector>(nullptr);
}

double GooCore::getTotalParticleMass(int idx) const
{
    double mass = particles_[idx]->mass;
    for (const auto& c : connectors_) {
        if (c->p1 == idx || c->p2 == idx)
            mass += 0.5 * c->mass;
    }
    return mass;
}

int GooCore::getNumRigidRods() const
{
    int nrods = 0;
    for (const auto& c : connectors_) {
        if (c->getType() == SimParameters::CT_RIGIDROD)
            nrods++;
    }
    return nrods;
}

std::tuple<Eigen::VectorXd, // q
           Eigen::VectorXd, // lambda
           Eigen::VectorXd> // v
GooCore::buildConfiguration() const
{
    Eigen::VectorXd q, lambda, v;
    int ndofs = 2 * particles_.size();
    q.resize(ndofs);
    v.resize(ndofs);

    for (int i = 0; i < (int)particles_.size(); i++) {
        q.segment<2>(2 * i) = particles_[i]->pos;
        v.segment<2>(2 * i) = particles_[i]->vel;
    }

    // TODO: Fill up the initial guessing of lambda
    // HINT: The dimension of lambda should depend on the number of rigid
    //       rods.
    lambda.setZero(getNumRigidRods());
    int i = 0;
    for (const auto& pc : connectors_) {
        if (pc->getType() == SimParameters::CT_RIGIDROD) {
            //std::shared_ptr<const RigidRod> pr = std::dynamic_pointer_cast<const RigidRod>(pc);
            lambda(i) = pc->loadLambda();
            i++;
        }
    }

    return std::make_tuple(q, lambda, v);
}

void GooCore::unbuildConfiguration(const Eigen::Ref<Eigen::VectorXd> q,
                                   const Eigen::Ref<Eigen::VectorXd> lambda,
                                   const Eigen::Ref<Eigen::VectorXd> v)
{
    int ndofs = q.size();
    assert(ndofs == int(2 * particles_.size()));
    int nrods = lambda.size();
    assert(nrods == getNumRigidRods());

    for (int i = 0; i < ndofs / 2; i++) {
        particles_[i]->pos = q.segment<2>(2 * i);
        particles_[i]->vel = v.segment<2>(2 * i);
    }

    // TODO: Save current lambda as the initial guessing of the next
    //       time step.
    int i = 0;
    for (auto& pc : connectors_) {
        if (pc->getType() == SimParameters::CT_RIGIDROD) {
            //std::shared_ptr<RigidRod> pr = std::dynamic_pointer_cast<RigidRod>(pc);
            pc->storeLambda(lambda, i);
            i++;
        }
    }
}

void GooCore::numericalIntegration(VectorXd& q, VectorXd& lambda, VectorXd& v)
{
    VectorXd F = createZeroForce();
    SparseMatrix<double> H;
    SparseMatrix<double> Minv;

    computeMassInverse(Minv);

    VectorXd oldq = q;

    switch (params_->constraintHandling) {
        case SimParameters::CH_PENALTY: {
            q += params_->timeStep * v;
            computeForceWithoutHessian(q, oldq, F);
            // TODO: Compute penalty force
            computePenaltyForce(q, F);
            v += params_->timeStep * Minv * F;
        } break;
        case SimParameters::CH_STEPPROJECT: {
            // TODO: Step and project
            // 1. "Step" with Velocity Verlet time integrator
            q += params_->timeStep * v;
            computeForceWithoutHessian(q, oldq, F);
            v += params_->timeStep * Minv * F;
            // 2. "Project" with projectOntoConstraints
            VectorXd unconstq = q;
            std::tie(q, lambda) = projectOntoConstraints(q, lambda, Minv);
            v += (q - unconstq) / params_->timeStep;
        } break;
        case SimParameters::CH_LAGRANGEMULT: {
            // TODO: Finish solveConstrainedLagrangian.
            (void)solveConstrainedLagrangian(q, lambda, v, Minv);
        } break;
    }
}

TripletArray GooCore::computeMassInverseCoeff() const
{
    int ndofs = 2 * int(particles_.size());
    std::vector<Eigen::Triplet<double>> Minvcoeffs;
    for (int i = 0; i < ndofs / 2; i++) {
        Minvcoeffs.emplace_back(2 * i, 2 * i, 1.0 / getTotalParticleMass(i));
        Minvcoeffs.emplace_back(2 * i + 1, 2 * i + 1, 1.0 / getTotalParticleMass(i));
    }
    return Minvcoeffs;
}

Eigen::VectorXd GooCore::createZeroForce() const
{
    Eigen::VectorXd f;
    int ndofs = 2 * particles_.size();
    f.setZero(ndofs);
    return f;
}

void GooCore::computeMassInverse(Eigen::SparseMatrix<double>& Minv)
{
    auto Minvcoeffs = computeMassInverseCoeff();
    int ndofs = 2 * int(particles_.size());
    Minv.resize(ndofs, ndofs);
    Minv.setFromTriplets(Minvcoeffs.begin(), Minvcoeffs.end());
}

void GooCore::computeForceWithoutHessian(const VectorXd& q,
                                         const VectorXd& qprev,
                                         Eigen::VectorXd& F) const
{
    F = createZeroForce();

    TripletArray Hcoeffs;
    if (params_->gravityEnabled)
        processGravityForce(F);
    if (params_->springsEnabled) {
        processSpringForce(q, F);
    }
    if (params_->dampingEnabled) {
        processDampingForce(q, qprev, F);
    }
    if (params_->floorEnabled) {
        processFloorForce(q, qprev, F);
    }
    if (params_->bendingEnabled) {
        processBendingForce(q, F);
    }
}

void GooCore::processGravityForce(Eigen::Ref<Eigen::VectorXd> F) const
{
    int nparticles = (int)particles_.size();
    for (int i = 0; i < nparticles; i++) {
        if (!particles_[i]->fixed) {
            F[2 * i + 1] += params_->gravityG * getTotalParticleMass(i);
        }
    }
}

TripletArray GooCore::processSpringForce(const Eigen::VectorXd& q,
                                         Eigen::Ref<Eigen::VectorXd> F) const
{
    TripletArray H;

    for (const auto& pc : connectors_) {
        pc->processSpringForce(q, F, H);
    }
    return H;
}

TripletArray GooCore::processDampingForce(const VectorXd& q,
                                          const VectorXd& qprev,
                                          Eigen::Ref<Eigen::VectorXd> F) const
{
    TripletArray H;

    for (const auto& pc : connectors_) {
        pc->processDampingForce(q, qprev, *params_, F, H);
    }
    return H;
}

TripletArray GooCore::processFloorForce(const VectorXd& q,
                                        const VectorXd& qprev,
                                        Eigen::Ref<Eigen::VectorXd> F) const
{
    int nparticles = particles_.size();

    double basestiffness = 10000;
    double basedrag = 1000.0;

    TripletArray H;

    for (int i = 0; i < nparticles; i++) {
        if (q[2 * i + 1] < -0.5 && !particles_[i]->fixed) {
            double vel = (q[2 * i + 1] - qprev[2 * i + 1]) / params_->timeStep;
            double dist = -0.5 - q[2 * i + 1];

            F[2 * i + 1] += basestiffness * dist - basedrag * dist * vel;

            double mag = basestiffness;
            mag += -0.5 * basedrag / params_->timeStep;
            mag += basedrag * qprev[2 * i + 1] / params_->timeStep;
            mag += -2.0 * basedrag * q[2 * i + 1] / params_->timeStep;
            H.emplace_back(2 * i + 1, 2 * i + 1, mag);
        }
    }
    return H;
}

void GooCore::processBendingForce(const Eigen::VectorXd& q,
                                  Eigen::Ref<Eigen::VectorXd> F) const
{
    // TODO: Bending force
    for (const auto & pb : bendingStencils_) {
        Vector3d x1, x2, x3;
        x1.segment<2>(0) = q.segment<2>(2 * pb.p1);
        x2.segment<2>(0) = q.segment<2>(2 * pb.p2);
        x3.segment<2>(0) = q.segment<2>(2 * pb.p3);
        x1(2) = x2(2) = x3(2) = 0;
        //std::cout << x1.transpose() << std::endl;
        //std::cout << x2.transpose() << std::endl;
        //std::cout << x3.transpose() << std::endl;
        Vector3d vec1 = x2 - x1, vec2 = x3 - x2, z_hat = Vector3d(0, 0, 1);
        double theta = 2 * atan2((vec1.cross(vec2)).dot(z_hat), vec1.norm() * vec2.norm() + vec1.dot(vec2));
        double k = pb.kb;
        Vector3d F1_3D = k * theta * (vec1.cross(z_hat)) / vec1.squaredNorm();
        Vector3d F3_3D = k * theta * (vec2.cross(z_hat)) / vec2.squaredNorm();
        Vector3d F2_3D = -F1_3D - F3_3D;
        Vector2d F1 = F1_3D.segment<2>(0), F2 = F2_3D.segment<2>(0), F3 = F3_3D.segment<2>(0);
        F.segment<2>(2 * pb.p1) += F1;
        F.segment<2>(2 * pb.p2) += F2;
        F.segment<2>(2 * pb.p3) += F3;
    }
}

void GooCore::computePenaltyForce(const Eigen::VectorXd& q,
                                  Eigen::Ref<Eigen::VectorXd> F) const
{
    // TODO: Implmentation of penalty force method
    for (const auto& pc : connectors_) {
        pc->computePenaltyForce(q, params_->penaltyStiffness, F);
    }
}

Eigen::VectorXd
GooCore::createZeroLagrangian() const
{
    int ndofs = 2 * particles_.size();
    int nrods = getNumRigidRods();
    Eigen::VectorXd f;
    f.setZero(ndofs + nrods);
    return f;
}

std::vector<std::shared_ptr<Connector>>
GooCore::listAllRigidRods()
{
    std::vector<std::shared_ptr<Connector>> ret;
    for (const auto& c : connectors_) {
        if (c->getType() == SimParameters::CT_RIGIDROD)
            ret.emplace_back(c);
    }
    return ret;
}

std::tuple<Eigen::VectorXd, // newq
           Eigen::VectorXd> // new lambda
GooCore::projectOntoConstraints(const VectorXd& q,
                                const Eigen::VectorXd& lambda,
                                const SparseMatrix<double>& Minv) const
{
    VectorXd guessq = q;
    VectorXd guesslambda = lambda;
    // TODO: Implment the "Project" step of the step and projection method
    VectorXd f = createZeroLagrangian(), dx = createZeroLagrangian();
    TripletArray Dfcoeffs = computeProjectionInfo(guessq, q, guesslambda, Minv, f);
    SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> solver;
    int numOfiter = 0;
    while (f.norm() > params_->NewtonTolerance && numOfiter < params_->NewtonMaxIters) {
        //std::cout << f.transpose() << std::endl;
        SparseMatrix<double> df(f.size(), f.size());
        df.setFromTriplets(Dfcoeffs.begin(), Dfcoeffs.end());
        //for (int i = 0; i < f.size();i++)
            //std::cout << df.row(i) << std::endl;
        //std::cout << std::endl;
        solver.compute(df);
        dx = solver.solve(-f);
        guessq += dx.segment(0, q.size());
        guesslambda += dx.segment(q.size(), lambda.size());
        Dfcoeffs = computeProjectionInfo(guessq, q, guesslambda, Minv, f);
        numOfiter++;
        if (numOfiter == params_->NewtonMaxIters)
            std::cout << "hit maxIter!" << std::endl;
    }

    return std::make_tuple(guessq, guesslambda);
}

TripletArray
GooCore::computeProjectionInfo(const Eigen::VectorXd& q,
                               const Eigen::VectorXd& qunconstrained,
                               const Eigen::VectorXd& lambda,
                               const Eigen::SparseMatrix<double>& Minv,
                               Eigen::Ref<Eigen::VectorXd> f) const
{
    std::vector<Eigen::Triplet<double>> Dfcoeffs;
    Dfcoeffs.reserve(q.size() + lambda.size() * 16);
    // TODO: Compute f and [df] for Step and Project Method.
    //f
    VectorXd g = createZeroConstraints();
    TripletArray dg = computeConstraintAndGradient(q, g);
    SparseMatrix<double> MatOfDg(g.size(), q.size());
    MatOfDg.setFromTriplets(dg.begin(), dg.end());
    SparseMatrix<double> MinvDgTrans = Minv * MatOfDg.transpose();
    VectorXd topNOfF = q - qunconstrained;
    for (int i = 0; i < lambda.size(); i++) {
        topNOfF += lambda(i) * MinvDgTrans.col(i);
    }
    f.segment(0, q.size()) = topNOfF;
    f.segment(q.size(), lambda.size()) = g;

    //df
    //for (int i = 0; i < q.size();i++)
        //std::cout << Minv.row(i) << std::endl;
    //std::cout << std::endl;
    for (int i = 0; i < dg.size(); i++) {
        Dfcoeffs.emplace_back(dg[i].row() + q.size(), dg[i].col(), dg[i].value());
        Dfcoeffs.emplace_back(dg[i].col(), dg[i].row() + q.size(), Minv.coeff(dg[i].col(), dg[i].col()) * dg[i].value());
        //std::cout << i << " " << Minv.coeff(dg[i].col(), dg[i].col()) << " " << dg[i].value() << std::endl;
    }
    for (int i = 0; i < q.size(); i++) {
        Dfcoeffs.emplace_back(i, i, 1);
    }
    for (int i = 0; i < lambda.size(); i++) {

    }
    //std::cout << "excuted!" << std::endl;
    for (const auto& pc : connectors_) {
        int i = 0;
        if (pc->getType() == SimParameters::CT_RIGIDROD) {
            Matrix2d I;
            I.setIdentity();
            Matrix2d localH = 2 * I;
            for (int j = 0; j < 2; j++) {
                for (int k = 0; k < 2; k++) {
                    int row_p1 = 2 * pc->p1 + j;
                    int col_p1 = 2 * pc->p1 + k;
                    int row_p2 = 2 * pc->p2 + j;
                    int col_p2 = 2 * pc->p2 + k;
                    Dfcoeffs.emplace_back(row_p1, col_p1, lambda(i) * Minv.coeff(row_p1, row_p1) * localH.coeff(j, k));
                    Dfcoeffs.emplace_back(row_p2, col_p2, lambda(i) * Minv.coeff(row_p2, row_p2) * localH.coeff(j, k));
                    Dfcoeffs.emplace_back(row_p1, col_p2, lambda(i) * Minv.coeff(row_p1, row_p1) * -localH.coeff(j, k));
                    Dfcoeffs.emplace_back(row_p2, col_p1, lambda(i) * Minv.coeff(row_p2, row_p2) * -localH.coeff(j, k));
                }
            }
            i++;
        }
    }
    return Dfcoeffs;
}

Eigen::VectorXd
GooCore::createZeroConstraints() const
{
    Eigen::VectorXd g;
    g.setZero(getNumRigidRods());
    return g;
}

TripletArray
GooCore::computeConstraintAndGradient(const VectorXd& q,
                                      Eigen::Ref<Eigen::VectorXd> g) const
{
    std::vector<Eigen::Triplet<double>> Dgcoeffs;
    Dgcoeffs.reserve(g.size() * 2);

    // TODO: Implement g and [dg]
    //int ndofs = q.size();
    int i = 0;
    for (const auto& pc : connectors_) {
        if (pc->getType() == SimParameters::CT_RIGIDROD) {
            Vector2d x1 = q.segment<2>(2 * pc->p1);
            Vector2d x2 = q.segment<2>(2 * pc->p2);
            std::shared_ptr<const RigidRod> pr = std::dynamic_pointer_cast<const RigidRod>(pc);
            g(i) = (x2 - x1).squaredNorm() - pr->length * pr->length;

            Vector2d localDg = 2 * (x2 - x1);
            for (int j = 0; j < 2; j++) {
                Dgcoeffs.emplace_back(i, 2 * pc->p1 + j, -localDg(j));
                //Dgcoeffs.emplace_back(2 * pc->p1 + j, ndofs + i, -localDg(j) / particles_[pc->p1]->mass);
                Dgcoeffs.emplace_back(i, 2 * pc->p2 + j, localDg(j));
                //Dgcoeffs.emplace_back(2 * pc->p2 + j, ndofs + i, localDg(j) / particles_[pc->p2]->mass);
            }

            i++;
        }
    }
    /*std::cout << q.transpose() << std::endl;
    std::cout << g.transpose() << std::endl;
    SparseMatrix<double> dg(g.size(), q.size());
    dg.setFromTriplets(Dgcoeffs.begin(), Dgcoeffs.end());
    for (int i = 0; i < g.size(); i++) {
        std::cout << dg.row(i) << std::endl;
    }*/
    return Dgcoeffs;
}

Eigen::SparseMatrix<double> // Df
GooCore::solveConstrainedLagrangian(Eigen::Ref<Eigen::VectorXd> q,
                                    Eigen::Ref<Eigen::VectorXd> lambda,
                                    Eigen::Ref<Eigen::VectorXd> v,
                                    const Eigen::SparseMatrix<double>& Minv) const
{
    Eigen::VectorXd oldq = q;
    Eigen::VectorXd lambdaguess = lambda;
    Eigen::VectorXd F;
    SparseMatrix<double> Df;
    int nrods = getNumRigidRods();
    // TODO: Implmentation of additional member functions
    double h = params_->timeStep;
    q += h * v;
    computeForceWithoutHessian(q, oldq, F);
    VectorXd g = createZeroConstraints();
    TripletArray Dgcoeffs_q = computeConstraintAndGradient(q, g);
    SparseMatrix<double> dg_q(lambda.size(), q.size());
    dg_q.setFromTriplets(Dgcoeffs_q.begin(), Dgcoeffs_q.end());
    VectorXd Q = q + h * v + h * h * Minv * F + h * h * Minv * dg_q.transpose() * lambdaguess;
    VectorXd f = createZeroConstraints();
    TripletArray Dgcoeffs_Q = computeConstraintAndGradient(Q, f);
    SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> solver;
    SparseMatrix<double> dg_Q(lambda.size(), q.size());
    int numOfiter = 0;
    while (f.norm() > params_->NewtonTolerance && numOfiter < params_->NewtonMaxIters) { 
        dg_Q.setFromTriplets(Dgcoeffs_Q.begin(), Dgcoeffs_Q.end());
        Df = h * h * dg_Q * Minv * dg_q.transpose();
        solver.compute(Df);
        lambdaguess += solver.solve(-f);
        Q = q + h * v + h * h * Minv * F + h * h * Minv * dg_q.transpose() * lambdaguess;
        Dgcoeffs_Q = computeConstraintAndGradient(Q, f);
        numOfiter++;
        if (numOfiter == params_->NewtonMaxIters)
            std::cout << "hit maxIter!" << std::endl;
    }
    v += Minv * (h * F + h * dg_q.transpose() * lambdaguess);
    return Df;
}

/****************************************************************
 * Objects removal.
 ****************************************************************/
double
GooCore::ptSegmentDist(const Vector2d& p,
                       const Vector2d& q1,
                       const Vector2d& q2) const
{
    double t = (p - q1).dot(q2 - q1) / (q2 - q1).dot(q2 - q1);
    t = std::max(0.0, std::min(t, 1.0));
    double linedistsq = (q1 + t * (q2 - q1) - p).squaredNorm();
    return sqrt(linedistsq);
}

void GooCore::detectSawedConnectors(std::set<int>& connectorsToDelete)
{
    for (int i = 0; i < (int)connectors_.size(); i++) {
        Vector2d pos1 = particles_[connectors_[i]->p1]->pos;
        Vector2d pos2 = particles_[connectors_[i]->p2]->pos;
        double maxx = std::max(pos1[0], pos2[0]);
        double minx = std::min(pos1[0], pos2[0]);
        double maxy = std::max(pos1[1], pos2[1]);
        double miny = std::min(pos1[1], pos2[1]);
        for (std::vector<Saw>::iterator saw = saws_.begin(); saw != saws_.end(); ++saw) {
            Vector2d sawpos = saw->pos;
            double sawr = saw->radius;

            if (sawpos[0] - sawr > maxx || sawpos[0] + sawr < minx || sawpos[1] - sawr > maxy || sawpos[1] + sawr < miny)
                continue;

            double sawspringdist = ptSegmentDist(sawpos, pos1, pos2);
            if (sawspringdist <= sawr) {
                connectorsToDelete.insert(i);
                break;
            }
        }
    }
}

void GooCore::detectSawedParticles(std::set<int>& particlesToDelete)
{
    for (int i = 0; i < (int)particles_.size(); i++) {
        Vector2d partpos = particles_[i]->pos;

        if (!(fabs(partpos[0]) < 2 && fabs(partpos[1]) < 2)) {
            particlesToDelete.insert(i);
            continue;
        }

        for (std::vector<Saw>::iterator it = saws_.begin(); it != saws_.end(); ++it) {
            Vector2d sawpos = it->pos;
            double sqdist = (sawpos - partpos).squaredNorm();
            if (sqdist < it->radius * it->radius) {
                particlesToDelete.insert(i);
                break;
            }
        }
    }
}

void GooCore::deleteSawedObjects()
{
    std::set<int> particlestodelete;
    std::set<int> connectorstodelete;
    detectSawedParticles(particlestodelete);
    detectSawedConnectors(connectorstodelete);

    decltype(particles_) newparticles;
    decltype(connectors_) newconnectors;
    std::vector<int> remainingparticlemap;

    std::set<int> bendingtodelete;
    std::vector<BendingStencil> newbending;
    std::vector<int> remainingbendingmap;

    if (!particlestodelete.empty()) {
        for (int i = 0; i < (int)connectors_.size(); i++) {
            if (particlestodelete.count(connectors_[i]->p1) || particlestodelete.count(connectors_[i]->p2))
                connectorstodelete.insert(i);
        }

        for (int i = 0; i < (int)particles_.size(); i++) {
            if (particlestodelete.count(i) == 0) {
                remainingparticlemap.emplace_back(newparticles.size());
                newparticles.emplace_back(particles_[i]);
            } else
                remainingparticlemap.emplace_back(-1);
        }
    }
    if (!connectorstodelete.empty()) {
        // TODO: Management of bending stencils, Part 1: bendingtodelete
        for (int i = 0; i < (int)connectors_.size(); i++) {
            if (connectorstodelete.count(i) == 0) {
                newconnectors.emplace_back(connectors_[i]);
            } else {
                if (!connectors_[i]->associatedBendingStencils.empty()) {
                    bendingtodelete.insert(connectors_[i]->associatedBendingStencils.begin(), connectors_[i]->associatedBendingStencils.end());
                }
                connectors_[i].reset();
            }
        }
    }
    // TODO: Management of bending stencils, Part 2: remainingbendingmap
    if (!bendingtodelete.empty()) {
        for (std::set<int>::iterator it = bendingtodelete.begin(); it != bendingtodelete.end(); ++it)
            std::cout << *it << " ";
        std::cout << std::endl;
        for (int i = 0; i < bendingStencils_.size(); i++) {
            if (bendingtodelete.count(i) == 0) {
                remainingbendingmap.emplace_back(newbending.size());
                newbending.emplace_back(bendingStencils_[i]);
            } else
                remainingbendingmap.emplace_back(-1);
        }
    }

    if (!connectorstodelete.empty() || !particlestodelete.empty()) {
        if (!connectorstodelete.empty())
            std::swap(connectors_, newconnectors);
    // TODO: Management of bending stencils, Part 3: update
    //       associatedBendingStencils
        if (!bendingtodelete.empty()) {
            std::swap(bendingStencils_, newbending);
            for (auto& c : connectors_) {
                std::set<int> newMap;
                for (std::set<int>::iterator it = c->associatedBendingStencils.begin(); it != c->associatedBendingStencils.end(); ++it) {
                    if (remainingbendingmap[*it]!=-1)
                        newMap.insert(remainingbendingmap[*it]);
                }
                c->associatedBendingStencils.swap(newMap);
            }
        }
        if (!particlestodelete.empty()) {
            std::swap(particles_, newparticles);
            for (auto& c : connectors_) {
                c->p1 = remainingparticlemap[c->p1];
                c->p2 = remainingparticlemap[c->p2];
            }
            for (auto& b : bendingStencils_) {
                b.p1 = remainingparticlemap[b.p1];
                b.p2 = remainingparticlemap[b.p2];
                b.p3 = remainingparticlemap[b.p3];
            }
            particle_offset_.clear();
            for (int i = 0; i < (int)particles_.size(); i++) {
                const auto& p = *particles_[i];
                particle_offset_[p.uid] = i;
            }
            // TODO: Management of bending stencils, Part 4: update
            //       p1, p2 and p3.
        }
    }
}

void GooCore::pruneOverstrainedSprings()
{
    int nsprings = connectors_.size();

    std::vector<int> toremove;
    for (int i = 0; i < nsprings; i++) {
        if (connectors_[i]->getType() != SimParameters::CT_SPRING)
            continue;

        Spring& s = *(Spring*)connectors_[i].get();
        if (s.canSnap) {
            Vector2d srcpos = particles_[s.p1]->pos;
            Vector2d dstpos = particles_[s.p2]->pos;
            double dist = (dstpos - srcpos).norm();

            double strain = (dist - s.restlen) / s.restlen;
            if (strain > params_->maxSpringStrain)
                toremove.emplace_back(i);
        }
    }

    for (std::vector<int>::reverse_iterator it = toremove.rbegin(); it != toremove.rend(); ++it) {
        assert(connectors_[*it]->associatedBendingStencils.empty());
        connectors_.erase(connectors_.begin() + *it);
    }
}

}
