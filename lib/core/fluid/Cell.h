#pragma once

#include <Eigen/Dense>
#include <unordered_set>

namespace fluid {

	enum CellType {
		EMPTY,
		FULL,
		SURFACE,
		BOUNDARY
	};

	class Cell {
	public:
		Cell(double size, Eigen::Vector2i index, Eigen::Vector2d pos, bool isBoundary);
		~Cell();

		void reset();
		void updateOldVel();

		double mSize, mPressure;
		Eigen::Vector2d mVelocity;
		Eigen::Vector2d mVelocity_old;
		bool mVelUIsKnown, mVelUIsKnown_old;
		bool mVelVIsKnown, mVelVIsKnown_old;
		Eigen::Vector2d mCenterPos;
		Eigen::Vector2i mCellIndex;
		CellType mType;
		Eigen::Vector4i mBoundaryCondision;//0 left, 1 bottom, 2 right, 3 top; 0 boundary, 1 water, 2 air;
		int mIndexInFulidCellArray;

		std::unordered_set<int> mMarkerParticles;
	};

	class MarkerParticle {
	public:
		MarkerParticle(Eigen::Vector2d position);
		~MarkerParticle();

		void reset();

		Eigen::Vector2d mPosition;
		Eigen::Vector2d mOriginalPos;
	};
}

