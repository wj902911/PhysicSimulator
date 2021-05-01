#include "Utils.h"
#include <igl/centroid.h>
#include <igl/per_face_normals.h>
#include <igl/doublearea.h>
#include <igl/massmatrix.h>

namespace cloth2 {
Eigen::MatrixXd centroids(Eigen::Ref<const Eigen::MatrixXd> V,
                          Eigen::Ref<const Eigen::MatrixXi> F)
{
    Eigen::MatrixXd C;

    C.resize(F.rows(), V.cols());
    for (int i = 0; i < F.rows(); i++) {
        Eigen::Vector3i f = F.row(i);
        C.row(i) = 1/3.0 * (V.row(f(0)) + V.row(f(1)) + V.row(f(2)));
    }
    return C;
}

Eigen::MatrixXd face_normals(Eigen::Ref<const Eigen::MatrixXd> V,
                             Eigen::Ref<const Eigen::MatrixXi> F)
{
    Eigen::MatrixXd N;
    igl::per_face_normals(V, F, N);
    return N;
}

Eigen::VectorXd face_areas(Eigen::Ref<const Eigen::MatrixXd> V,
                           Eigen::Ref<const Eigen::MatrixXi> F)
{
    Eigen::VectorXd A;
    igl::doublearea(V, F, A);
    A *= 0.5;
    return A;
}

Eigen::SparseMatrix<double> inv_mass_matrix(Eigen::Ref<const Eigen::MatrixXd> V,
                                            Eigen::Ref<const Eigen::MatrixXi> F)
{
    Eigen::SparseMatrix<double> M;
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
    return M.cwiseInverse();
}

}
