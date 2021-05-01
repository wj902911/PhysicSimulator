#include "ClothCore.h"
#include "SimParameters.h"
#include "SceneObjects.h"
#include "Utils.h"
#include "../shared/VectorMath.h"
#include "CollisionDetection.h"
#include "debug_IntermediateRecord.h"
#include <Eigen/SVD>
#include <iostream>
#include <algorithm>
#include <limits>
#include <map>
#include <set>
#include <unordered_set>
#include <Eigen/Dense>

namespace cloth2 {

IntermediateRecord::IntermediateRecord()
{
}

IntermediateRecord::~IntermediateRecord()
{
}

ClothCore::ClothCore()
{
    params_ = std::make_shared<SimParameters>();
    clickedVertex_ = -1;
}

ClothCore::~ClothCore()
{
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd>
ClothCore::getCurrentMesh() const
{
    int totverts = Q_.rows();
    int totfaces = F_.rows();
    for (auto &it : bodyInstances_) {
        totverts += it->getTemplate().getVerts().rows();
        totfaces += it->getTemplate().getFaces().rows();
    }

    Eigen::MatrixXd renderQ(totverts, 3);
    int clothverts = Q_.rows();
    int vertidx = 0;
    for (int i = 0; i < clothverts; i++) {
        renderQ.row(vertidx++) = Q_.row(i);
    }
    std::vector<int> vertoffsets;
    for (auto& it : bodyInstances_) {
        vertoffsets.push_back(vertidx);
        int nbodyverts = it->getTemplate().getVerts().rows();
        for (int i = 0; i < nbodyverts; i++) {
            Eigen::Vector3d templatePt = it->getTemplate().getVerts().row(i).transpose();
            Eigen::Vector3d worldPt = VectorMath::rotationMatrix(it->theta) * templatePt + it->c;
            renderQ.row(vertidx++) = worldPt;
        }
    }
    Eigen::MatrixXi renderF(totfaces, 3);
    int faceidx = 0;
    int clothfaces = F_.rows();
    for (int i = 0; i < clothfaces; i++) {
        renderF.row(faceidx++) = F_.row(i);
    }
    int bodyidx = 0;
    for (auto& it : bodyInstances_) {
        int nbodyfaces = it->getTemplate().getFaces().rows();
        for (int i=0; i<nbodyfaces; i++) {
            for (int j=0; j<3; j++) {
                renderF(faceidx, j) = it->getTemplate().getFaces()(i, j) + vertoffsets[bodyidx];
            }
            faceidx++;
        }
        bodyidx++;
    }

    Eigen::MatrixXd C;
    return std::make_tuple(renderQ, renderF, C);
}

void ClothCore::attachMesh(Eigen::Ref<Eigen::MatrixX3d> V,
                           Eigen::Ref<Eigen::MatrixX3i> F,
                           double scale)
{
    F_ = F;
    origF_ = F;
    origQ_ = V * scale;
    // rotate mesh
    Eigen::Matrix3d rot;
    rot << -1.0, 0.0, 0.0,
        0.0, 0.0, 1.0,
        0.0, -1.0, 0.0;
    origQ_ = origQ_ * rot.transpose();
    // translate
    origQ_.col(1).setConstant(1.0);
    initSimulation();
}

void ClothCore::attachBody(Eigen::Ref<Eigen::MatrixX3d> V,
                           Eigen::Ref<Eigen::MatrixX3i> F,
                           double scale)
{
    body_ = std::make_shared<RigidBodyTemplate>(V, F, scale);
    launchBody_ = std::make_shared<RigidBodyTemplate>(V, F, 0.25 * scale);
    initSimulation();
}

void ClothCore::initSimulation()
{
    Q_ = origQ_;
    F_ = origF_;
    Qdot_.resize(Q_.rows(), 3);
    Qdot_.setZero();
    pinnedVerts_.clear();
    bodyInstances_.clear();
    params_->shouldDeployBevo = false;
    bevoDeployed_ = false;


    int nverts = Q_.rows();
    int topleft = -1;
    int topright = -1;
    double topleftdist = -std::numeric_limits<double>::infinity();
    double toprightdist = -std::numeric_limits<double>::infinity();
    Eigen::Vector3d tr(1, 1, 0);
    Eigen::Vector3d tl(-1, 1, 0);
    for (int i = 0; i < nverts; i++) {
        double disttr = tr.dot(Q_.row(i));
        if (disttr > toprightdist) {
            toprightdist = disttr;
            topright = i;
        }
        double disttl = tl.dot(Q_.row(i));
        if (disttl > topleftdist) {
            topleftdist = disttl;
            topleft = i;
        }
    }
    pinnedVerts_.push_back(topleft);
    pinnedVerts_.push_back(topright);

    Minv_ = inv_mass_matrix(origQ_, F_);
    computeHinges();
    // TODO: Implmentation of additional member functions
}

bool ClothCore::simulateOneStep()
{
    int nverts = Q_.rows();

    // update rigid bodies

    std::vector<Eigen::Vector3d> oldthetas;
    for (int bodyidx = 0; bodyidx < (int)bodyInstances_.size(); bodyidx++) {
        RigidBodyInstance& body = *bodyInstances_[bodyidx];
        body.oldc = body.c;
        body.c += params_->dt * body.cvel;
        Eigen::Matrix3d Rhw = VectorMath::rotationMatrix(params_->dt * body.w);
        Eigen::Matrix3d Rtheta = VectorMath::rotationMatrix(body.theta);

        Eigen::Vector3d oldtheta = body.theta;
        body.theta = VectorMath::axisAngle(Rtheta * Rhw);
        if (body.theta.dot(oldtheta) < 0 && oldtheta.norm() > M_PI / 2.0) {
            double oldnorm = oldtheta.norm();
            oldtheta = (oldnorm - 2.0 * M_PI) * oldtheta / oldnorm;
        }

        oldthetas.push_back(oldtheta);
        body.oldtheta = oldtheta;
    }

    // rigid body precession

    for (int bodyidx = 0; bodyidx < (int)bodyInstances_.size(); bodyidx++) {
        RigidBodyInstance& body = *bodyInstances_[bodyidx];
        Eigen::Matrix3d Mi = body.getTemplate().getInertiaTensor();

        Eigen::Vector3d newwguess(body.w);

        int iter = 0;
        for (iter = 0; iter < 20; iter++) {
            Eigen::Vector3d term1 = (-VectorMath::TMatrix(-params_->dt * newwguess).inverse() * VectorMath::TMatrix(oldthetas[bodyidx])).transpose() * Mi * body.density * newwguess;
            Eigen::Vector3d term2 = (VectorMath::TMatrix(params_->dt * body.w).inverse() * VectorMath::TMatrix(oldthetas[bodyidx])).transpose() * Mi * body.density * body.w;
            Eigen::Vector3d fval = term1 + term2;
            if (fval.norm() / body.density / Mi.trace() <= 1e-6)
                break;

            Eigen::Matrix3d Df = (-VectorMath::TMatrix(-params_->dt * newwguess).inverse() * VectorMath::TMatrix(body.theta)).transpose() * Mi * body.density;

            Eigen::Vector3d deltaw = Df.inverse() * (-fval);
            newwguess += deltaw;
        }
        // std::cout << "Converged in " << iter << " Newton iterations" << std::endl;
        body.w = newwguess;
    }

    Eigen::MatrixXd oldQ = Q_;
    Q_ += params_->dt * Qdot_;

    if (params_->enableIntermediateValueRecording) {
        record_oldQ = oldQ;
        records.resize(params_->constraintIters);
    }
    // apply constraints
    for (int i = 0; i < params_->constraintIters; i++) {
        applyConstraints();

        // collision detection

        if (params_->collisionsEnabled) {
            // TODO: Implement findCollisions
            auto col_info = findCollisions(oldQ, Q_, F_, bodyInstances_);

            if (params_->enableIntermediateValueRecording) {
                records[i].Q = Q_;
                records[i].col_info = col_info;
            }

            // TODO: applyClothVertexVsBodyFaceConstraints
            applyClothVertexVsBodyFaceConstraints(col_info);
            if (params_->enableIntermediateValueRecording) {
                records[i].Q_after_cvbf = Q_;
            }
            // TODO: applyClothFaceVsBodyVertexConstraints
            applyClothFaceVsBodyVertexConstraints(col_info);
            if (params_->enableIntermediateValueRecording) {
                records[i].Q_after_cfbv = Q_;
            }
        }
    }

    Qdot_ = (Q_ - oldQ) / params_->dt;

    if (params_->gravityEnabled) {
        Qdot_.col(1) = Qdot_.col(1).array() + params_->dt * params_->gravityG;
    }
    Qdot_ += params_->dt * Minv_ * computeJetForces();

    return false;
}

void ClothCore::hold(int vertex_id, const Eigen::Vector3d& position)
{
    // std::cerr << "hold vertex " << vertex_id << " at " << position.transpose() << std::endl;
    clickedVertex_ = vertex_id;
    curPos_ = position;
}

void ClothCore::updateHold(const Eigen::Vector3d& position)
{
    if (clickedVertex_ < 0)
        return;
    // std::cerr << "update hold vertex " << " to " << position.transpose() << std::endl;
    curPos_ = position;
}

void ClothCore::releaseHold()
{
    clickedVertex_ = -1;
}


void ClothCore::cutFaces(const std::vector<int>& faces)
{
    if (faces.empty())
        return;
    std::unordered_set<int> todelete;

    for (auto it : faces)
        if (it >= 0 && it < F_.rows())
            todelete.insert(it);

    Eigen::MatrixXi newF(F_.rows() - todelete.size(), F_.cols());
    int new_index = 0;
    for (int i = 0; i < F_.rows(); i++) {
        if (todelete.count(i) == 0) {
            newF.row(new_index) = F_.row(i);
            new_index++;
        }
    }
    std::swap(F_, newF);
    Minv_ = inv_mass_matrix(origQ_, F_);
    computeHinges();
}

void ClothCore::computeHinges()
{
    std::map<std::pair<int, int>, std::vector<int>> edgemap;
    int nfaces = F_.rows();
    for (int i = 0; i < nfaces; i++) {
        for (int j = 0; j < 3; j++) {
            int nextj = (j + 1) % 3;
            int v1 = F_(i, j);
            int v2 = F_(i, nextj);
            if (v1 > v2)
                std::swap(v1, v2);
            edgemap[std::pair<int, int>(v1, v2)].push_back(i);
        }
    }

    int nhinges = 0;
    for (auto it : edgemap) {
        if (it.second.size() == 2)
            nhinges++;
    }
    H_.resize(nhinges, 4);
    int idx = 0;
    for (auto it : edgemap) {
        if (it.second.size() != 2)
            continue;
        std::set<int> hingeverts;
        for (int j = 0; j < 3; j++) {
            hingeverts.insert(F_(it.second[0], j));
            hingeverts.insert(F_(it.second[1], j));
        }
        int colidx = 0;
        for (auto v : hingeverts) {
            H_(idx, colidx) = v;
            colidx++;
        }
        idx++;
    }
}

void ClothCore::applyConstraints()
{
    if (params_->pinEnabled) {
        applyPinConstraints();
    }

    if (params_->stretchEnabled) {
        applyStretchConstraints();
    }

    if (params_->bendingEnabled) {
        applyBendingConstraints();
    }

    if (clickedVertex_ != -1 && params_->dragEnabled) {
        applyDragConstraints(clickedVertex_);
    }
}

void ClothCore::applyPinConstraints()
{
    for (int idx : pinnedVerts_) {
        Q_.row(idx) = params_->pinWeight * origQ_.row(idx) + (1.0 - params_->pinWeight) * Q_.row(idx);
    }
}

std::tuple<Eigen::Vector3d, // oldcentroid
           Eigen::Vector3d, // newcentroid
           Eigen::Matrix3d> // The solution
ClothCore::solveOrthogonalProcrustes(Eigen::VectorXi primitive) const
{
    int nelem = primitive.size();
    Eigen::Vector3d oldcentroid(0, 0, 0);
    Eigen::Vector3d newcentroid(0, 0, 0);
    for (int j = 0; j < nelem; j++) {
        oldcentroid += origQ_.row(primitive(j));
        newcentroid += Q_.row(primitive(j));
    }
    oldcentroid /= nelem;
    newcentroid /= nelem;
    Eigen::Matrix<double, 3, -1> A(3, nelem);
    Eigen::Matrix<double, 3, -1> B(3, nelem);
    for (int j = 0; j < nelem; j++) {
        Eigen::Vector3d oldvec = origQ_.row(primitive(j)).transpose() - oldcentroid;
        A.col(j) = oldvec;
        Eigen::Vector3d newvec = Q_.row(primitive(j)).transpose() - newcentroid;
        B.col(j) = newvec;
    }
    Eigen::Matrix3d M = B * A.transpose();
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d R = svd.matrixU() * svd.matrixV().transpose();
    return std::make_tuple(oldcentroid, newcentroid, R);
}

void ClothCore::applyStretchConstraints()
{
    int nfaces = F_.rows();
    for (int i = 0; i < nfaces; i++) {
        auto tup = solveOrthogonalProcrustes(F_.row(i));
        auto oldcentroid = std::get<0>(tup);
        auto newcentroid = std::get<1>(tup);
        auto R = std::get<2>(tup);
        for (int j = 0; j < 3; j++) {
            Eigen::Vector3d oldvec = origQ_.row(F_(i, j)).transpose() - oldcentroid;
            Eigen::Vector3d newpt = R * oldvec + newcentroid;
            Q_.row(F_(i, j)) = params_->stretchWeight * newpt.transpose() + (1 - params_->stretchWeight) * Q_.row(F_(i, j));
        }
    }
}

void ClothCore::applyBendingConstraints()
{
    int nhinges = H_.rows();
    for (int i = 0; i < nhinges; i++) {
        auto tup = solveOrthogonalProcrustes(H_.row(i));
        auto oldcentroid = std::get<0>(tup);
        auto newcentroid = std::get<1>(tup);
        auto R = std::get<2>(tup);
        for (int j = 0; j < 4; j++) {
            Eigen::Vector3d oldvec = origQ_.row(H_(i, j)).transpose() - oldcentroid;
            Eigen::Vector3d newpt = R * oldvec + newcentroid;
            Q_.row(H_(i, j)) = params_->bendingWeight * newpt.transpose() + (1 - params_->bendingWeight) * Q_.row(H_(i, j));
        }
    }
}

void ClothCore::applyDragConstraints(int vertex)
{
    if (clickedVertex_ < 0 || clickedVertex_ >= Q_.rows())
        return;

    Q_.row(clickedVertex_) = params_->dragWeight * curPos_.transpose() + (1.0 - params_->dragWeight) * Q_.row(clickedVertex_);
}

void ClothCore::applyClothVertexVsBodyFaceConstraints(const ColInfo& col_info)
{
    if (!params_->enableCVBF)
        return ;
    const auto& vertCollisions = col_info.vertCollisions;
    // TODO: Implmentation Body Vertex vs. Cloth Face Contact Response
    for (int i = 0; i < vertCollisions.size(); i++) {
        if (vertCollisions[i].empty())
            continue;
        Eigen::Vector3d q0e = Q_.row(i);
        const auto& face_set = vertCollisions[i];
        for (const auto& face_pair : face_set) {
            int body_index = face_pair.first;
            int face_index = face_pair.second;
            Eigen::Vector3i face = bodyInstances_[body_index]->getTemplate().getFaces().row(face_index);
            Eigen::Vector3d q1s, q1e, q2s, q2e, q3s, q3e;
            std::tie(q1s, q1e) = bodyInstances_[body_index]->getVertPos(face(0));
            std::tie(q2s, q2e) = bodyInstances_[body_index]->getVertPos(face(1));
            std::tie(q3s, q3e) = bodyInstances_[body_index]->getVertPos(face(2));
            Eigen::Vector3d normal_end = (q2e - q1e).cross(q3e - q1e).normalized();
            double projOnNorm = (q0e - q1e).dot(normal_end);
            if (projOnNorm < 0) {
                Q_.row(i) = q0e - (projOnNorm - 1e-6) * normal_end;
            }
        }
    }
}

void ClothCore::applyClothFaceVsBodyVertexConstraints(const ColInfo& col_info)
{
    if (!params_->enableCFBV)
        return ;
    const auto& faceCollisions = col_info.faceCollisions;
    // TODO: Implmentation Cloth Face vs. Body Vertex Contact Response
    for (int i = 0; i < faceCollisions.size(); i++) {
        if (faceCollisions[i].empty())
            continue;
        Eigen::Vector3i face = F_.row(i);
        Eigen::Vector3d q1e = Q_.row(face(0));
        Eigen::Vector3d q2e = Q_.row(face(1));
        Eigen::Vector3d q3e = Q_.row(face(2));
        const auto& vert_set = faceCollisions[i];
        for (const auto& vert_pair : vert_set) {
            Eigen::Vector3d q0e, q0s;
            std::tie(q0s, q0e) = bodyInstances_[vert_pair.first]->getVertPos(vert_pair.second);
            Eigen::Vector3d normal_end = (q2e - q1e).cross(q3e - q1e).normalized();
            double projOnNorm = (q0e - q1e).dot(normal_end);
            if (projOnNorm < 0) {
                Q_.row(face(0)) = q1e + (projOnNorm - 1e-6) * normal_end;
                Q_.row(face(1)) = q2e + (projOnNorm - 1e-6) * normal_end;
                Q_.row(face(2)) = q3e + (projOnNorm - 1e-6) * normal_end;
            }
        }
    }
}


uint32_t
ClothCore::addJetSource(const Eigen::Vector3d& xc, // Source
                        const Eigen::Vector3d& dir, // Direction
                        double jetStrength,
                        double jetAtten)
{
    auto ret = jets_id_;
    jets_.emplace_back(std::make_shared<JetSource>(xc, dir, jetStrength, jetAtten, ret));
    jets_id_++;

    return ret;
}

std::shared_ptr<JetSource>
ClothCore::queryJetSource(uint32_t id)
{
    for (const auto j : jets_)
        if (j->id() == id)
            return j;

    return std::shared_ptr<JetSource>(nullptr);
}

void ClothCore::removeJetSource(uint32_t id)
{
    auto func = [id](const std::shared_ptr<JetSource>& j) {
        return j->id() == id;
    };
    jets_.erase(std::remove_if(jets_.begin(), jets_.end(), func),
                jets_.end());
}

void ClothCore::clearJetSource()
{
    jets_.clear();
}

Eigen::MatrixXd
ClothCore::computeJetForces() const
{
    Eigen::MatrixXd forces;
    forces.setZero(Q_.rows(), Q_.cols());

    Eigen::MatrixXd C, N;
    Eigen::VectorXd A;
    // TODO: Implement centroids to compute the centroids from Q_
    C = centroids(Q_, F_);
    N = face_normals(Q_, F_);
    A = face_areas(Q_, F_);

    for (const auto& pj : jets_) {
        // TODO: Implement this function
        pj->blow(Q_, F_, A, N, C, forces);
    }
    return forces;
}

void ClothCore::deployBevo()
{
    if (bevoDeployed_)
        return;

    Eigen::Vector3d center(0, 0, 0);
    Eigen::Vector3d orientation(0, 0, 0);
    Eigen::Vector3d cvel(0, 0, 0);
    Eigen::Vector3d omega(0, 0, 0);

    bodyInstances_.push_back(std::make_shared<RigidBodyInstance>(body_, center, orientation, cvel, omega, 1.0));

    bevoDeployed_ = true;
}

void ClothCore::launchBevo(Eigen::Vector3d src,
                           Eigen::Vector3d dir,
                           Eigen::Vector3d orientation,
                           Eigen::Vector3d omega)
{
    bodyInstances_.push_back(std::make_shared<RigidBodyInstance>(launchBody_,
                                                                 src + 0.1*dir,
                                                                 orientation,
                                                                 params_->launchSpeed * dir,
                                                                 params_->launchSpinningSpeed * omega,
                                                                 1.0));
}

}
