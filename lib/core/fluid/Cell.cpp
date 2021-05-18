#include "Cell.h"

namespace fluid {

	Cell::Cell(double size, Eigen::Vector2i index, Eigen::Vector2d pos, bool isBoundary)
	{
		mSize = size;
		mCellIndex = index;
		mPressure = 0;
		mVelocity = mVelocity_old = Eigen::Vector2d::Zero();
		mVelUIsKnown = false;
		mVelUIsKnown_old = true;
		mVelVIsKnown = false;
		mVelVIsKnown_old = true;
		mCenterPos = pos;
		if (isBoundary) {
			mType = CellType::BOUNDARY;
			mBoundaryCondision = Eigen::Vector4i::Zero();
		}
		else {
			mType = CellType::EMPTY;
			mBoundaryCondision = Eigen::Vector4i::Ones();
		}
		mIndexInFulidCellArray = -1;
	}

	Cell::~Cell()
	{
	}

	void Cell::reset()
	{
		mPressure = 0;
		mVelocity = Eigen::Vector2d::Zero();
		if (mType != CellType::BOUNDARY) {
			mBoundaryCondision = Eigen::Vector4i::Ones();
			mType = CellType::EMPTY;
		}
		mMarkerParticles.clear();
	}

	void Cell::updateOldVel()
	{
		mVelocity_old = mVelocity;
		mVelUIsKnown_old = mVelUIsKnown;
		mVelVIsKnown_old = mVelVIsKnown;
	}

	MarkerParticle::MarkerParticle(Eigen::Vector2d position)
	{
		mPosition = position;
		mOriginalPos = position;
	}

	MarkerParticle::~MarkerParticle()
	{
	}

	void MarkerParticle::reset()
	{
		mPosition = mOriginalPos;
	}
}


