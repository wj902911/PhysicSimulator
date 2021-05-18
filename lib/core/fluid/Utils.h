#pragma once

#include <eigen/Eigen/Core>

namespace fluid {

	double bilinearInterpolation(Eigen::Vector2d pos, double x1, double x2, double y1, double y2, double value1, double value2, double value3, double value4);

}