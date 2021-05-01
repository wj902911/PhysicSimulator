#ifndef PSIM_CORE_CLOTH2_DEBUG_INTERMEDIATE_RECORD_H
#define PSIM_CORE_CLOTH2_DEBUG_INTERMEDIATE_RECORD_H

namespace cloth2 {

#include <Eigen/Core>
#include "CollisionDetection.h"

struct IntermediateRecord {
    IntermediateRecord();
    ~IntermediateRecord();

    Eigen::MatrixXd Q;
    ColInfo col_info;
    Eigen::MatrixXd Q_after_cvbf;
    Eigen::MatrixXd Q_after_cfbv;
};

}

#endif
