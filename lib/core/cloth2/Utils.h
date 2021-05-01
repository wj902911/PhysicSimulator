#ifndef PSIM_CORE_CLOTH2_UTILS_H
#define PSIM_CORE_CLOTH2_UTILS_H

#include <Eigen/Core>
#include <Eigen/SparseCore>


namespace cloth2 {
Eigen::MatrixXd centroids(Eigen::Ref<const Eigen::MatrixXd> V,
                          Eigen::Ref<const Eigen::MatrixXi> F);

Eigen::MatrixXd face_normals(Eigen::Ref<const Eigen::MatrixXd> V,
                             Eigen::Ref<const Eigen::MatrixXi> F);

Eigen::VectorXd face_areas(Eigen::Ref<const Eigen::MatrixXd> V,
                           Eigen::Ref<const Eigen::MatrixXi> F);

Eigen::SparseMatrix<double> inv_mass_matrix(Eigen::Ref<const Eigen::MatrixXd> V,
                                            Eigen::Ref<const Eigen::MatrixXi> F);

}

#endif
