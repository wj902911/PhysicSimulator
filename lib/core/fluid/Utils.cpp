#include "Utils.h"

namespace fluid {
	double bilinearInterpolation(Eigen::Vector2d pos, double x1, double x2, double y1, double y2, double value1, double value2, double value3, double value4)
	{
		double x = pos(0), y = pos(1);
		double A = (x2 - x1) * (y2 - y1);
		double A1 = (x2 - x) * (y2 - y);
		double A2 = (x - x1) * (y2 - y);
		double A3 = (x - x1) * (y - y1);
		double A4 = (x2 - x) * (y - y1);
		return (value1 * A1 + value2 * A2 + value3 * A3 + value4 * A4) / A;
	}
}