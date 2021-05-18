#include "FluidCore.h"
#include "Utils.h"
#include <iostream>
#include <fstream>

namespace fluid {

	FluidCore::FluidCore()
	{
		mParams_= std::make_shared<SimParameters>();
		mCellSize = 0.1;
		mSubTimeStep = mParams_->timeStep;

		initModel();
	}

	FluidCore::~FluidCore()
	{
	}

	void FluidCore::initSimulation()
	{
		for (int i = 0; i < mCells_.size(); i++) {
			for (int j = 0; j < mCells_.size(); j++) {
				mCells_[i][j]->reset();
			}
		}

		for (int i = 0; i < mMarkerParticles_.size(); i++)
			mMarkerParticles_[i]->reset();

		mUActiveCells_.clear();
		mVActiveCells_.clear();
		mFluidCells_.clear();

		assignParticle();
		flagCell();
		//initPressure();
	}

	bool FluidCore::simulateOneStep()
	{
		assignParticle();
		flagCell();
		buildActiveCellArray();
		//buildConfiguration();
		determineTimeStep();
		advection();
		applyBodyForce();
		initPressure();
		project();
		setBoundaryVel();
		//outputVelField();
		updateOldVelocity();
		updateParticlePos();

		return false;
	}

	std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd> 
		FluidCore::getCurrentMesh() const
	{
		double baseradius = 0.025;
		int numcirclewedges = 20;

		std::vector<Eigen::Vector3d> verts;
		std::vector<Eigen::Vector3d> vertexColors;
		std::vector<Eigen::Vector3i> faces;

		int idx = 0;

		double eps = 1e-4;

		for (int i = 0; i < 24; i++) {
			vertexColors.emplace_back(0.1, 0.1, 0.1);
		}

		int cellrows = mCells_.size();
		int cellcols = mCells_.size();
		double thickness=0.05;
		Eigen::Vector3d bl, bll, blb, br, brr, brb, tr, trr, trt, tl, tll, tlt;

		bl = Eigen::Vector3d(mCells_[0][0]->mCenterPos(0) + mCellSize * 0.5, mCells_[0][0]->mCenterPos(1) + mCellSize * 0.5, eps);
		bll = Eigen::Vector3d(bl(0) - thickness, bl(1), bl(2));
		blb = Eigen::Vector3d(bl(0), bl(1) - thickness, bl(2));

		br = Eigen::Vector3d(mCells_[cellcols - 1][0]->mCenterPos(0) - mCellSize * 0.5, mCells_[cellcols - 1][0]->mCenterPos(1) + mCellSize * 0.5, eps);
		brr = Eigen::Vector3d(br(0) + thickness, br(1), br(2));
		brb = Eigen::Vector3d(br(0), br(1) - thickness, br(2));

		tr = Eigen::Vector3d(mCells_[cellcols - 1][cellrows - 1]->mCenterPos(0) - mCellSize * 0.5, mCells_[cellcols - 1][cellrows - 1]->mCenterPos(1) - mCellSize * 0.5, eps);
		trr = Eigen::Vector3d(tr(0) + thickness, tr(1), tr(2));
		trt = Eigen::Vector3d(tr(0), tr(1) + thickness, tr(2));

		tl = Eigen::Vector3d(mCells_[0][cellrows - 1]->mCenterPos(0) + mCellSize * 0.5, mCells_[0][cellrows - 1]->mCenterPos(1) - mCellSize * 0.5, eps);
		tll = Eigen::Vector3d(tl(0) - thickness, tl(1), tl(2));
		tlt = Eigen::Vector3d(tl(0), tl(1) + thickness, tl(2));

		verts.emplace_back(tl);
		verts.emplace_back(tll);
		verts.emplace_back(bll);
		faces.emplace_back(idx, idx + 1, idx + 2);

		verts.emplace_back(bll);
		verts.emplace_back(bl);
		verts.emplace_back(tl);
		faces.emplace_back(idx + 3, idx + 4, idx + 5);
		idx += 6;

		verts.emplace_back(bl);
		verts.emplace_back(blb);
		verts.emplace_back(br);
		faces.emplace_back(idx, idx + 1, idx + 2);

		verts.emplace_back(blb);
		verts.emplace_back(brb);
		verts.emplace_back(br);
		faces.emplace_back(idx + 3, idx + 4, idx + 5);
		idx += 6;

		verts.emplace_back(br);
		verts.emplace_back(brr);
		verts.emplace_back(trr);
		faces.emplace_back(idx, idx + 1, idx + 2);

		verts.emplace_back(trr);
		verts.emplace_back(tr);
		verts.emplace_back(br);
		faces.emplace_back(idx + 3, idx + 4, idx + 5);
		idx += 6;

		verts.emplace_back(tlt);
		verts.emplace_back(tl);
		verts.emplace_back(tr);
		faces.emplace_back(idx, idx + 1, idx + 2);

		verts.emplace_back(tr);
		verts.emplace_back(trt);
		verts.emplace_back(tlt);
		faces.emplace_back(idx + 3, idx + 4, idx + 5);
		idx += 6;

		int nparticles = mMarkerParticles_.size();

		for (int i = 0; i < nparticles; i++) {
			double radius = baseradius;

			Eigen::Vector3d color(0, 1, 0);

			for (int j = 0; j < numcirclewedges + 2; j++) {
				vertexColors.emplace_back(color);
			}

			verts.emplace_back(mMarkerParticles_[i]->mPosition[0], mMarkerParticles_[i]->mPosition[1], 0);

			const double PI = 3.1415926535898;
			for (int j = 0; j <= numcirclewedges; j++) {
				verts.emplace_back(mMarkerParticles_[i]->mPosition[0] + radius * cos(2 * PI * j / numcirclewedges),
					mMarkerParticles_[i]->mPosition[1] + radius * sin(2 * PI * j / numcirclewedges),
					0);
			}

			for (int j = 0; j <= numcirclewedges; j++) {
				faces.emplace_back(idx, idx + j + 1, idx + 1 + ((j + 1) % (numcirclewedges + 1)));
			}

			idx += numcirclewedges + 2;
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

	std::tuple<Eigen::VectorXd,
		Eigen::VectorXd,
		Eigen::VectorXd>
		FluidCore::buildConfiguration() const
	{
		Eigen::VectorXd q, p, v;
		int nCells = mCells_[0].size() * mCells_.size();
		int nParticles = mMarkerParticles_.size();
		q.resize(2 * nParticles);
		p.resize(nCells);
		v.resize(2 * nCells);

		int n = 0;
		for (int i = 0; i < mCells_.size(); i++) {
			for (int j = 0; j < mCells_[i].size(); j++) {
				if (mCells_[i][j]->mType != CellType::BOUNDARY && mCells_[i][j]->mType != CellType::EMPTY) {
					v.segment<2>(2 * n) = mCells_[i][j]->mVelocity_old;
					p(n) = mCells_[i][j]->mPressure;
					n++;
				}
			}
		}

		for (int i = 0; i < nParticles; i++) {
			q.segment<2>(2 * i) = mMarkerParticles_[i]->mPosition;
		}

		return std::make_tuple(q, p, v);
	}

	Eigen::Vector2i FluidCore::mapPosToCell(Eigen::Vector2d pos, double cellSize)
	{
		Eigen::Vector2i index = floor((pos / cellSize).array()).cast<int>();
		return index;
	}

	Eigen::Vector2i FluidCore::mapPosToInterpolatCell_v(Eigen::Vector2d pos, double cellSize)
	{
		int i = floor(pos(0) / cellSize - 0.5), j = floor(pos(1) / cellSize);
		Eigen::Vector2i index(i, j);
		return index;
	}

	Eigen::Vector2i FluidCore::mapPosToInterpolatCell_u(Eigen::Vector2d pos, double cellSize)
	{
		int i = floor(pos(0) / cellSize), j = floor(pos(1) / cellSize - 0.5);
		Eigen::Vector2i index(i, j);
		return index;
	}

	void FluidCore::initModel()
	{
		int fluidCellrows = 9, fluidCellcols = 9;
		for (int i = 0; i < fluidCellcols; i++) {
			for (int j = 0; j < fluidCellrows; j++) {
				Eigen::MatrixX2d ppos(4, 2);
				ppos << (i + 1) * mCellSize + 0.025, (j + 1) * mCellSize + 0.025,
					    (i + 1) * mCellSize + 0.075, (j + 1) * mCellSize + 0.025,
					    (i + 1) * mCellSize + 0.025, (j + 1) * mCellSize + 0.075,
					    (i + 1) * mCellSize + 0.075, (j + 1) * mCellSize + 0.075;
				        //(i + 1) * mCellSize + 0.05, (j + 1) * mCellSize + 0.05;
				for (int k = 0; k < ppos.rows(); k++) {
					mMarkerParticles_.emplace_back(std::make_shared<MarkerParticle>(MarkerParticle(ppos.row(k))));
				}
			}
		}
		
		int cellrows = 20, cellcols = 20;
		for (int i = 0; i < cellrows; i++) {
			std::vector<std::shared_ptr<Cell>> cellcol;
			for (int j = 0; j < cellcols; j++) {
				Eigen::Vector2i index(i, j);
				Eigen::Vector2d cpos((i + 0.5) * mCellSize, (j + 0.5) * mCellSize);
				bool isBoundary;
				if (i == 0 || i == cellcols - 1 || j == 0 || j == cellrows - 1) {
					isBoundary = true;
				}
				else
					isBoundary = false;
				cellcol.emplace_back(std::make_shared<Cell>(Cell(mCellSize, index, cpos, isBoundary)));
			}
			mCells_.emplace_back(cellcol);
		}
	}

	void FluidCore::assignParticle()
	{
		for (int i = 0; i < mFluidCells_.size(); i++) {
			mFluidCells_[i]->mMarkerParticles.clear();
		}
		for (int i = 0; i < mMarkerParticles_.size(); i++) {
			Eigen::Vector2i index = mapPosToCell(mMarkerParticles_[i]->mPosition, mCellSize);
			mCells_[index(0)][index(1)]->mMarkerParticles.insert(i);
		}
	}

	void FluidCore::flagCell()
	{
		for (int i = 0; i < mCells_.size(); i++) {
			for (int j = 0; j < mCells_[i].size(); j++) {
				if (mCells_[i][j]->mType != CellType::BOUNDARY) {
					if (mCells_[i][j]->mMarkerParticles.size() != 0) {
						bool adjacentEmpty = false;
						for (int k = 0; k < 2; k++) {
							if (mCells_[i - 1 + 2 * k][j]->mType != CellType::BOUNDARY) {
								if (mCells_[i - 1 + 2 * k][j]->mMarkerParticles.size() == 0) {
									adjacentEmpty = true;
									mCells_[i][j]->mBoundaryCondision(2 * k) = 2;
								}
								else
									mCells_[i][j]->mBoundaryCondision(2 * k) = 1;
							}
							else {
								mCells_[i][j]->mBoundaryCondision(2 * k) = 0;
							}
							if (mCells_[i][j - 1 + 2 * k]->mType != CellType::BOUNDARY) {
								if (mCells_[i][j - 1 + 2 * k]->mMarkerParticles.size() == 0) {
									adjacentEmpty = true;
									mCells_[i][j]->mBoundaryCondision(1 + 2 * k) = 2;
								}
								else
									mCells_[i][j]->mBoundaryCondision(1 + 2 * k) = 1;
							}
							else {
								mCells_[i][j]->mBoundaryCondision(1 + 2 * k) = 0;
							}
						}
						if(adjacentEmpty)
							mCells_[i][j]->mType = CellType::SURFACE;
						else
							mCells_[i][j]->mType = CellType::FULL;
					}
					else {
						mCells_[i][j]->mType = CellType::EMPTY;
						mCells_[i][j]->mIndexInFulidCellArray = -1;
					}
				}
			}
		}
	}

	void FluidCore::initPressure()
	{
#if 0
		for (int i = 0; i < mCells_.size(); i++) {
			int n = 0;
			for (int j = mCells_[i].size(); j > 0; j--) {
				if (mCells_[i][j - 1]->mType != CellType::BOUNDARY) {
					if (mCells_[i][j - 1]->mMarkerParticles.size() > 0) {
						mCells_[i][j - 1]->mPressure = mParams_->airPressure + (n + 1) * mParams_->density * mParams_->gravity * mCellSize;
						n++;
					}
					else {
						mCells_[i][j - 1]->mPressure = mParams_->airPressure;
						n = 0;
					}
				}
			}
		}
#endif
		for (int i = 0; i < mCells_.size(); i++) {
			for (int j =  0; j < mCells_[i].size(); j++) {
				mCells_[i][j]->mPressure = 0;
			}
		}
	}

	void FluidCore::buildActiveCellArray()
	{
		mUActiveCells_.clear();
		mVActiveCells_.clear();
		mFluidCells_.clear();

		for (int i = 0; i < mCells_.size(); i++) {
			for (int j = 0; j < mCells_[i].size(); j++) {
				if (mCells_[i][j]->mType != CellType::BOUNDARY) {
					if (mCells_[i][j]->mType == CellType::EMPTY) {
						if (mCells_[i - 1][j]->mMarkerParticles.size() > 0) {
							mUActiveCells_.emplace_back(mCells_[i][j]);
							mCells_[i][j]->mVelUIsKnown = true;
						}
						else
							mCells_[i][j]->mVelUIsKnown = false;

						if (mCells_[i][j - 1]->mMarkerParticles.size() > 0) {
							mVActiveCells_.emplace_back(mCells_[i][j]);
							mCells_[i][j]->mVelVIsKnown = true;
						}
						else
							mCells_[i][j]->mVelVIsKnown = false;
					}
					else {
						mCells_[i][j]->mIndexInFulidCellArray = mFluidCells_.size();
						mFluidCells_.emplace_back(mCells_[i][j]);

						mCells_[i][j]->mVelUIsKnown = true;
						mCells_[i][j]->mVelVIsKnown = true;

						if (mCells_[i - 1][j]->mType != CellType::BOUNDARY) {
							mUActiveCells_.emplace_back(mCells_[i][j]);
						}
						else
							mCells_[i][j]->mVelocity(0) = 0;

						if (mCells_[i][j - 1]->mType != CellType::BOUNDARY) {
							mVActiveCells_.emplace_back(mCells_[i][j]);
						}
						else
							mCells_[i][j]->mVelocity(1) = 0;
					}
				}
#if 1
				else {
					if (i == mCells_.size() - 1) {
						if (mCells_[i - 1][j]->mMarkerParticles.size() > 0) {
							mCells_[i][j]->mVelocity(0) = 0;
							mCells_[i][j]->mVelUIsKnown = true;
						}
						else
							mCells_[i][j]->mVelUIsKnown = false;
					}

					if (j == mCells_[i].size() - 1) {
						if (mCells_[i][j - 1]->mMarkerParticles.size() > 0) {
							mCells_[i][j]->mVelocity(1) = 0;
							mCells_[i][j]->mVelVIsKnown = true;
						}
						else
							mCells_[i][j]->mVelVIsKnown = false;
					}
				}
#endif
			}
		}
	}

	void FluidCore::determineTimeStep()
	{
		Eigen::VectorXd velMag = mParticleVels_.rowwise().norm();
		double maxSpeed;
		if (velMag.size() == 0)
			maxSpeed = 0;
		else
			maxSpeed= velMag.maxCoeff();
		double maxDt = mCellSize / (maxSpeed + 10e-30);
		mSubTimeStep = mParams_->timeStep < maxDt ? mParams_->timeStep : maxDt;
	}

	//advection
	void FluidCore::advection()
	{
		double h = mSubTimeStep;
		for (int i = 0; i < mUActiveCells_.size(); i++) {
			//Eigen::Vector2i index = mUActiveCells_[i]->mCellIndex;
			double dx = mUActiveCells_[i]->mSize;
			Eigen::Vector2d xg = mUActiveCells_[i]->mCenterPos;
			xg(0) -= 0.5 * dx;
			Eigen::Vector2d xmid = xg - 0.5 * h * getCellVelAtU(mUActiveCells_[i]);
			Eigen::Vector2d xp = xg - h * getPosVel(xmid, dx);
			mUActiveCells_[i]->mVelocity(0) = getPosVel_u(xp, dx);
		}
		for (int i = 0; i < mVActiveCells_.size(); i++) {
			double dx = mVActiveCells_[i]->mSize;
			Eigen::Vector2d xg = mVActiveCells_[i]->mCenterPos;
			xg(1) -= 0.5 * dx;
			Eigen::Vector2d xmid = xg - 0.5 * h * getCellVelAtV(mVActiveCells_[i]);
			Eigen::Vector2d xp = xg - h * getPosVel(xmid, dx);
			mVActiveCells_[i]->mVelocity(1) = getPosVel_v(xp, dx);
		}
	}

	void FluidCore::applyBodyForce()
	{
		double h = mSubTimeStep;
		double g = mParams_->gravity;
		for (int i = 0; i < mVActiveCells_.size(); i++) {
			mVActiveCells_[i]->mVelocity(1) += h * g;
		}
	}

	void FluidCore::project()
	{
		Eigen::SparseMatrix<double> A(mFluidCells_.size(), mFluidCells_.size());
		TripletArray coeffA = getPressureLaplacianCoeff();
		A.setFromTriplets(coeffA.begin(), coeffA.end());
		Eigen::VectorXd b = getRHS();
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
		solver.compute(A);
		if (solver.info() != Eigen::Success) {
			std::cerr << "compute(A) failed!" << std::endl;
			return;
		}
		Eigen::VectorXd pressure = mParams_->density * solver.solve(b) / mSubTimeStep;
		if (solver.info() != Eigen::Success) {
			std::cerr << "solve(b) failed!" << std::endl;
			return;
		}

		for (int i = 0; i < mFluidCells_.size(); i++) {
			mFluidCells_[i]->mPressure = pressure[i];
		}

		for (int i = 0; i < mUActiveCells_.size(); i++) {
			Eigen::Vector2i index = mUActiveCells_[i]->mCellIndex;
			mUActiveCells_[i]->mVelocity(0) -= mSubTimeStep / mParams_->density / mCellSize * (mUActiveCells_[i]->mPressure - mCells_[index(0) - 1][index(1)]->mPressure);
		}
		for (int i = 0; i < mVActiveCells_.size(); i++) {
			Eigen::Vector2i index = mVActiveCells_[i]->mCellIndex;
			mVActiveCells_[i]->mVelocity(1) -= mSubTimeStep / mParams_->density / mCellSize * (mVActiveCells_[i]->mPressure - mCells_[index(0)][index(1) - 1]->mPressure);
		}
	}

	void FluidCore::setBoundaryVel()
	{
		int sign = 1;
		for (int i = 0; i < mCells_.size(); i++) {
			for (int j = 0; j < mCells_[i].size(); j++) {
				if (mCells_[i][j]->mType == CellType::BOUNDARY) {
					if (i == mCells_.size() - 1) {
						if (mCells_[i - 1][j]->mVelVIsKnown) {
							mCells_[i][j]->mVelocity(1) = sign * mCells_[i - 1][j]->mVelocity(1);
							mCells_[i][j]->mVelVIsKnown = true;
						}
						else
							mCells_[i][j]->mVelVIsKnown = false;
					}

					if (i == 0 ) {
						if (mCells_[i + 1][j]->mVelVIsKnown) {
							mCells_[i][j]->mVelocity(1) = sign * mCells_[i + 1][j]->mVelocity(1);
							mCells_[i][j]->mVelVIsKnown = true;
						}
						else
							mCells_[i][j]->mVelVIsKnown = false;
					}

					if (j == mCells_[i].size() - 1) {
						if (mCells_[i][j - 1]->mVelUIsKnown) {
							mCells_[i][j]->mVelocity(0) = sign * mCells_[i][j - 1]->mVelocity(0);
							mCells_[i][j]->mVelUIsKnown = true;
						}
						else
							mCells_[i][j]->mVelUIsKnown = false;
					}

					if (j == 0) {
						if (mCells_[i][j + 1]->mVelUIsKnown) {
							mCells_[i][j]->mVelocity(0) = sign * mCells_[i][j + 1]->mVelocity(0);
							mCells_[i][j]->mVelUIsKnown = true;
						}
						else
							mCells_[i][j]->mVelUIsKnown = false;
					}
				}
#if 0
				if (mCells_[i][j]->mType != CellType::BOUNDARY) {
					if (mCells_[i][j]->mMarkerParticles.size()>0) {
						if (mCells_[i - 1][j]->mType == CellType::BOUNDARY) {
							mCells_[i - 1][j]->mVelocity(1) = -mCells_[i][j]->mVelocity(1);
							mCells_[i - 1][j]->mVelVIsKnown = true;
						}
						else
							mCells_[i - 1][j]->mVelVIsKnown = false;

						if (mCells_[i][j - 1]->mType == CellType::BOUNDARY) {
							mCells_[i][j - 1]->mVelocity(1) = -mCells_[i][j]->mVelocity(0);
							mCells_[i][j - 1]->mVelVIsKnown = true;
						}
						else
							mCells_[i][j - 1]->mVelVIsKnown = false;
					}
				}
				else {
					if (i != 0 && mCells_[i - 1][j]->mMarkerParticles.size() > 0) {
						mCells_[i][j]->mVelocity(0) = 0;
						mCells_[i][j]->mVelocity(1) = -mCells_[i - 1][j]->mVelocity(1);
						mCells_[i][j]->mVelUIsKnown = true;
						mCells_[i][j]->mVelVIsKnown = true;
					}
					else {
						mCells_[i][j]->mVelUIsKnown = false;
						mCells_[i][j]->mVelVIsKnown = false;
					}

					if (j != 0 && mCells_[i][j - 1]->mMarkerParticles.size() > 0) {
						mCells_[i][j]->mVelocity(1) = 0;
						mCells_[i][j]->mVelocity(0) = -mCells_[i][j - 1]->mVelocity(0);
						mCells_[i][j]->mVelVIsKnown = true;
						mCells_[i][j]->mVelUIsKnown = true;
					}
					else {
						mCells_[i][j]->mVelUIsKnown = false;
						mCells_[i][j]->mVelVIsKnown = false;
					}
				}
#endif
			}
		}
	}

	void FluidCore::updateOldVelocity()
	{
		for (int i = 0; i < mCells_.size(); i++) {
			for (int j = 0; j < mCells_[i].size(); j++) {
				mCells_[i][j]->updateOldVel();
			}
		}
	}

	void FluidCore::updateParticlePos()
	{
		mParticleVels_.resize(mMarkerParticles_.size(), 2);
		double h = mSubTimeStep;
		for (int i = 0; i < mMarkerParticles_.size(); i++) {
			Eigen::Vector2d xg = mMarkerParticles_[i]->mPosition;
			Eigen::Vector2d xmid = xg + 0.5 * h * getPosVel(xg, mCellSize);
			Eigen::Vector2d velMid = getPosVel(xmid, mCellSize);
			Eigen::Vector2d xc = xg + h * velMid;
			mMarkerParticles_[i]->mPosition = xc;
			mParticleVels_.row(i) = velMid;
		}
	}

	TripletArray FluidCore::getPressureLaplacianCoeff() const
	{
		TripletArray Coeff;
		//int m = 0, n = 0;
		for (int i = 0; i < mFluidCells_.size(); i++) {
			for (int j = 0; j < 4; j++) {
				if (mFluidCells_[i]->mBoundaryCondision(j) != 0) {
					Coeff.emplace_back(i, i, -1);
					if (mFluidCells_[i]->mBoundaryCondision(j) == 1) {
						Eigen::Vector2i index = mFluidCells_[i]->mCellIndex;
						if (j == 0 || j == 2)
							Coeff.emplace_back(i, mCells_[index(0) - 1 + j][index(1)]->mIndexInFulidCellArray, 1);
						else
							Coeff.emplace_back(i, mCells_[index(0)][index(1) - 2 + j]->mIndexInFulidCellArray, 1);
					}
				}
			}
		}
		return Coeff;
	}

	Eigen::VectorXd FluidCore::getRHS() const
	{
		Eigen::VectorXd rhs(mFluidCells_.size());
		for (int i = 0; i < mFluidCells_.size(); i++) {
			double u1, u2, v1, v2;
			Eigen::Vector2i index = mFluidCells_[i]->mCellIndex;

			if (mFluidCells_[i]->mBoundaryCondision(0) != 0)
				u1 = mFluidCells_[i]->mVelocity(0);
			else
				u1 = 0;

			if (mFluidCells_[i]->mBoundaryCondision(2) != 0)
				u2 = mCells_[index(0) + 1][index(1)]->mVelocity(0);
			else
				u2 = 0;

			if (mFluidCells_[i]->mBoundaryCondision(1) != 0)
				v1 = mFluidCells_[i]->mVelocity(1);
			else
				v1 = 0;

			if (mFluidCells_[i]->mBoundaryCondision(3) != 0)
				v2 = mCells_[index(0)][index(1) + 1]->mVelocity(1);
			else
				v2 = 0;

			rhs(i) = u2 - u1 + v2 - v1;
		}
		rhs *= mCellSize;
		return rhs;
	}

	Eigen::Vector2d FluidCore::getCellVelAtU(std::shared_ptr<Cell> cell)
	{
		Eigen::Vector2i index = cell->mCellIndex;
		double u = 0.0;
		if(cell->mVelUIsKnown_old)
			u = cell->mVelocity_old(0);
		else 
			u = extrapolateUVelOnGrid(index);

		double v = 0.0;
		if (!getVVelAtU(index, v)) {
			int n = 0;
			if (getVVelAtU(Eigen::Vector2i(index(0) - 1, index(1) + 0), v))
				n++;
			if (getVVelAtU(Eigen::Vector2i(index(0) + 0, index(1) - 1), v))
				n++;
			if (getVVelAtU(Eigen::Vector2i(index(0) + 1, index(1) + 0), v))
				n++;
			if (getVVelAtU(Eigen::Vector2i(index(0) + 0, index(1) + 1), v))
				n++;
			if (n != 0)
				v /= n;
			else
				std::cerr << "Cell:" << index << "Can not find known v at u for adv!" << std::endl;
		}
		Eigen::Vector2d vel(u, v);
		return vel;
	}

	Eigen::Vector2d FluidCore::getCellVelAtV(std::shared_ptr<Cell> cell)
	{
		Eigen::Vector2i index = cell->mCellIndex;
		double v = 0.0;
		if (cell->mVelVIsKnown_old)
			v = cell->mVelocity_old(1);
		else
			v = extrapolateVVelOnGrid(index);

		double u = 0.0;
		if (!getUVelAtV(index, u)) {
			int n = 0;
			if (getUVelAtV(Eigen::Vector2i(index(0) - 1, index(1) + 0), u))
				n++;
			if (getUVelAtV(Eigen::Vector2i(index(0) + 0, index(1) - 1), u))
				n++;
			if (getUVelAtV(Eigen::Vector2i(index(0) + 1, index(1) + 0), u))
				n++;
			if (getUVelAtV(Eigen::Vector2i(index(0) + 0, index(1) + 1), u))
				n++;
			if (n != 0)
				u /= n;
			else
				std::cerr << "Cell:" << index << "Can not find known u at v for adv!" << std::endl;
		}
		Eigen::Vector2d vel(u, v);
		return vel;
	}

	double FluidCore::extrapolateUVelOnGrid(Eigen::Vector2i index)
	{
		int n = 0;
		double valueSum = 0;

		if (index(0) != 0 && mCells_[index(0) - 1][index(1) + 0]->mVelUIsKnown_old) {
			n++;
			valueSum += mCells_[index(0) - 1][index(1) + 0]->mVelocity_old(0);
		}
		if (index(1) != 0 && mCells_[index(0) + 0][index(1) - 1]->mVelUIsKnown_old) {
			n++;
			valueSum += mCells_[index(0) + 0][index(1) - 1]->mVelocity_old(0);
		}
		if (index(0) != mCells_.size() - 1 && mCells_[index(0) + 1][index(1) + 0]->mVelUIsKnown_old) {
			n++;
			valueSum += mCells_[index(0) + 1][index(1) + 0]->mVelocity_old(0);
		}
		if (index(1) != mCells_[index(0)].size() - 1 && mCells_[index(0) + 0][index(1) + 1]->mVelUIsKnown_old) {
			n++;
			valueSum += mCells_[index(0) + 0][index(1) + 1]->mVelocity_old(0);
		}

		if (n == 0) {
			//std::cerr << "Can not find known u for adv!" << std::endl;
			return 0.0;
		}
		else
			return valueSum / n;
	}

	double FluidCore::extrapolateVVelOnGrid(Eigen::Vector2i index)
	{
		int n = 0;
		double valueSum = 0;

		if (index(0) != 0 && mCells_[index(0) - 1][index(1) + 0]->mVelVIsKnown_old) {
			n++;
			valueSum += mCells_[index(0) - 1][index(1) + 0]->mVelocity_old(1);
		}
		if (index(1) != 0 && mCells_[index(0) + 0][index(1) - 1]->mVelVIsKnown_old) {
			n++;
			valueSum += mCells_[index(0) + 0][index(1) - 1]->mVelocity_old(1);
		}
		if (index(0) != mCells_.size() - 1 && mCells_[index(0) + 1][index(1) + 0]->mVelVIsKnown_old) {
			n++;
			valueSum += mCells_[index(0) + 1][index(1) + 0]->mVelocity_old(1);
		}
		if (index(1) != mCells_[index(0)].size() - 1 && mCells_[index(0) + 0][index(1) + 1]->mVelVIsKnown_old) {
			n++;
			valueSum += mCells_[index(0) + 0][index(1) + 1]->mVelocity_old(1);
		}

		if (n == 0) {
			//std::cerr << "Can not find known v for adv!" << std::endl;
			return 0.0;
		}
		else
			return valueSum / n;
	}

	bool FluidCore::getVVelAtU(Eigen::Vector2i index, double& v)
	{
		int n = 0;
		double valueSum = 0;

		if (index(0) != 0 && mCells_[index(0) - 1][index(1) + 0]->mVelVIsKnown_old) {
			n++;
			valueSum += mCells_[index(0) - 1][index(1) + 0]->mVelocity_old(1);
		}
		if (mCells_[index(0) + 0][index(1) + 0]->mVelVIsKnown_old) {
			n++;
			valueSum += mCells_[index(0) + 0][index(1) + 0]->mVelocity_old(1);
		}
		if (index(0) != 0 && index(1) != mCells_[index(0)].size() - 1 && mCells_[index(0) - 1][index(1) + 1]->mVelVIsKnown_old) {
			n++;
			valueSum += mCells_[index(0) - 1][index(1) + 1]->mVelocity_old(1);
		}
		if (index(1) != mCells_[index(0)].size() - 1 && mCells_[index(0) + 0][index(1) + 1]->mVelVIsKnown_old) {
			n++;
			valueSum += mCells_[index(0) + 0][index(1) + 1]->mVelocity_old(1);
		}

		if (n == 0) {
			return false;
		}
		else {
			v += valueSum / n;
			return true;
		}
	}

	bool FluidCore::getUVelAtV(Eigen::Vector2i index, double& u)
	{
		int n = 0;
		double valueSum = 0;

		if (index(1) != 0 && mCells_[index(0)][index(1) - 1]->mVelUIsKnown_old) {
			n++;
			valueSum += mCells_[index(0)][index(1) - 1]->mVelocity_old(0);
		}
		if (index(0) != mCells_.size() - 1 && index(1) != 0 && mCells_[index(0) + 1][index(1) - 1]->mVelUIsKnown_old) {
			n++;
			valueSum += mCells_[index(0) + 1][index(1) - 1]->mVelocity_old(0);
		}
		if (index(0) != mCells_.size() - 1 && mCells_[index(0) + 1][index(1)]->mVelUIsKnown_old) {
			n++;
			valueSum += mCells_[index(0) + 1][index(1)]->mVelocity_old(0);
		}
		if (mCells_[index(0)][index(1)]->mVelUIsKnown_old) {
			n++;
			valueSum += mCells_[index(0)][index(1)]->mVelocity_old(0);
		}

		if (n == 0) {
			return false;
		}
		else {
			u += valueSum / n;
			return true;
		}
	}

	Eigen::Vector2d FluidCore::getPosVel(Eigen::Vector2d pos, double cellSize)
	{
		//u
		Eigen::Vector2i index_Corner = mapPosToInterpolatCell_u(pos, cellSize);
		Eigen::Vector2d centerPos = mCells_[index_Corner(0)][index_Corner(1)]->mCenterPos;
		centerPos(1) += 0.5 * cellSize;
		double u = getVelCompU(index_Corner, centerPos, pos, cellSize);
		//v
		index_Corner = mapPosToInterpolatCell_v(pos, cellSize);
		centerPos = mCells_[index_Corner(0)][index_Corner(1)]->mCenterPos;
		centerPos(0) += 0.5 * cellSize;
		double v = getVelCompV(index_Corner, centerPos, pos, cellSize);
		return Eigen::Vector2d(u, v);
	}

	double FluidCore::getPosVel_u(Eigen::Vector2d pos, double cellSize)
	{
		Eigen::Vector2i index_Corner = mapPosToInterpolatCell_u(pos, cellSize);
		Eigen::Vector2d centerPos = mCells_[index_Corner(0)][index_Corner(1)]->mCenterPos;
		centerPos(1) += 0.5 * cellSize;
		return getVelCompU(index_Corner, centerPos, pos, cellSize);
	}

	double FluidCore::getPosVel_v(Eigen::Vector2d pos, double cellSize)
	{
		Eigen::Vector2i index_Corner = mapPosToInterpolatCell_v(pos, cellSize);
		Eigen::Vector2d centerPos = mCells_[index_Corner(0)][index_Corner(1)]->mCenterPos;
		centerPos(0) += 0.5 * cellSize;
		return getVelCompV(index_Corner, centerPos, pos, cellSize);
	}

	double FluidCore::getVelCompU(Eigen::Vector2i index_Corner, Eigen::Vector2d centerPos, Eigen::Vector2d pointPos, double cellSize)
	{
		double x1 = centerPos(0) - 0.5 * cellSize;
		double x2 = centerPos(0) + 0.5 * cellSize;
		double y1 = centerPos(1) - 0.5 * cellSize;
		double y2 = centerPos(1) + 0.5 * cellSize;

		double value1, value2, value3, value4;

		if (mCells_[index_Corner(0) + 0][index_Corner(1) + 0]->mVelUIsKnown_old)
			value1 = mCells_[index_Corner(0) + 0][index_Corner(1) + 0]->mVelocity_old(0);
		else
			value1 = extrapolateUVelOnGrid(Eigen::Vector2i(index_Corner(0) + 0, index_Corner(1) + 0));

		if (mCells_[index_Corner(0) + 1][index_Corner(1) + 0]->mVelUIsKnown_old)
			value2 = mCells_[index_Corner(0) + 1][index_Corner(1) + 0]->mVelocity_old(0);
		else
			value2 = extrapolateUVelOnGrid(Eigen::Vector2i(index_Corner(0) + 1, index_Corner(1) + 0));

		if (mCells_[index_Corner(0) + 1][index_Corner(1) + 1]->mVelUIsKnown_old)
			value3 = mCells_[index_Corner(0) + 1][index_Corner(1) + 1]->mVelocity_old(0);
		else
			value3 = extrapolateUVelOnGrid(Eigen::Vector2i(index_Corner(0) + 1, index_Corner(1) + 1));

		if (mCells_[index_Corner(0) + 0][index_Corner(1) + 1]->mVelUIsKnown_old)
			value4 = mCells_[index_Corner(0) + 0][index_Corner(1) + 1]->mVelocity_old(0);
		else
			value4 = extrapolateUVelOnGrid(Eigen::Vector2i(index_Corner(0) + 0, index_Corner(1) + 1));

		return bilinearInterpolation(pointPos, x1, x2, y1, y2, value1, value2, value3, value4);
	}

	double FluidCore::getVelCompV(Eigen::Vector2i index_Corner, Eigen::Vector2d centerPos, Eigen::Vector2d pointPos, double cellSize)
	{
		double x1 = centerPos(0) - 0.5 * cellSize;
		double x2 = centerPos(0) + 0.5 * cellSize;
		double y1 = centerPos(1) - 0.5 * cellSize;
		double y2 = centerPos(1) + 0.5 * cellSize;

		double value1, value2, value3, value4;

		if (mCells_[index_Corner(0) + 0][index_Corner(1) + 0]->mVelVIsKnown_old)
			value1 = mCells_[index_Corner(0) + 0][index_Corner(1) + 0]->mVelocity_old(1);
		else
			value1 = extrapolateVVelOnGrid(Eigen::Vector2i(index_Corner(0) + 0, index_Corner(1) + 0));

		if (mCells_[index_Corner(0) + 1][index_Corner(1) + 0]->mVelVIsKnown_old)
			value2 = mCells_[index_Corner(0) + 1][index_Corner(1) + 0]->mVelocity_old(1);
		else
			value2 = extrapolateVVelOnGrid(Eigen::Vector2i(index_Corner(0) + 1, index_Corner(1) + 0));

		if (mCells_[index_Corner(0) + 1][index_Corner(1) + 1]->mVelVIsKnown_old)
			value3 = mCells_[index_Corner(0) + 1][index_Corner(1) + 1]->mVelocity_old(1);
		else
			value3 = extrapolateVVelOnGrid(Eigen::Vector2i(index_Corner(0) + 1, index_Corner(1) + 1));

		if (mCells_[index_Corner(0) + 0][index_Corner(1) + 1]->mVelVIsKnown_old)
			value4 = mCells_[index_Corner(0) + 0][index_Corner(1) + 1]->mVelocity_old(1);
		else
			value4 = extrapolateVVelOnGrid(Eigen::Vector2i(index_Corner(0) + 0, index_Corner(1) + 1));

		return bilinearInterpolation(pointPos, x1, x2, y1, y2, value1, value2, value3, value4);
	}

	void FluidCore::outputVelField()
	{
		std::ofstream fout1("VelField_U.txt"),fout2("VelField_V.txt");
		for (int i = 0; i < mCells_.size(); i++) {
			for (int j = 0; j < mCells_[i].size(); j++) {
				if (mCells_[j][mCells_[i].size() - 1 - i]->mVelUIsKnown)
					fout1 << mCells_[j][mCells_[i].size() - 1 - i]->mVelocity(0) << " ";
				else
					fout1 << "unknow ";

				if (mCells_[j][mCells_[i].size() - 1 - i]->mVelVIsKnown)
					fout2 << mCells_[j][mCells_[i].size() - 1 - i]->mVelocity(1) << " ";
				else
					fout2 << "unknow ";
			}
			fout1 << std::endl;
			fout2 << std::endl;
		}
		fout1.close();
		fout2.close();
	}
}
