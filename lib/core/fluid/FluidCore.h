#pragma once

#include "../PhysicsCore.h"
#include "Cell.h"
#include "SimParameters.h"
#include <Eigen/StdVector>
#include <Eigen/Sparse>


namespace fluid {

    using TripletArray = std::vector<Eigen::Triplet<double>>;

    class FluidCore : public PhysicsCore {
    public:
        FluidCore();
        ~FluidCore();

        void initSimulation() override;

        bool simulateOneStep() override;

        std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd>
            getCurrentMesh() const override;

        std::shared_ptr<SimParameters> getPointerToSimParameters()
        {
            return mParams_;
        }

        std::tuple<Eigen::VectorXd, // marker position (q)
                   Eigen::VectorXd, // pressure (p)
                   Eigen::VectorXd> // velocity (v)
        buildConfiguration() const;

        Eigen::Vector2i mapPosToCell(Eigen::Vector2d pos, double cellSize);
        Eigen::Vector2i mapPosToInterpolatCell_v(Eigen::Vector2d pos, double cellSize);
        Eigen::Vector2i mapPosToInterpolatCell_u(Eigen::Vector2d pos, double cellSize);

        void initModel();
        
        void assignParticle();

        void flagCell();

        void initPressure();

        void buildActiveCellArray();

        //void buildExtrapolateCellArray();
        void determineTimeStep();

        void advection();

        void applyBodyForce();

        void project();

        //void updatePressure();

        void setBoundaryVel();

        void updateOldVelocity();

        void updateParticlePos();

        TripletArray getPressureLaplacianCoeff() const;
        Eigen::VectorXd getRHS() const;

        Eigen::Vector2d getCellVelAtU(std::shared_ptr<Cell> cell);
        Eigen::Vector2d getCellVelAtV(std::shared_ptr<Cell> cell);

        double extrapolateUVelOnGrid(Eigen::Vector2i index);
        double extrapolateVVelOnGrid(Eigen::Vector2i index);
        bool getVVelAtU(Eigen::Vector2i index, double& v);
        bool getUVelAtV(Eigen::Vector2i index, double& u);

        Eigen::Vector2d getPosVel(Eigen::Vector2d pos, double cellSize);
        double getPosVel_u(Eigen::Vector2d pos, double cellSize);
        double getPosVel_v(Eigen::Vector2d pos, double cellSize);
        double getVelCompU(Eigen::Vector2i index_Corner, Eigen::Vector2d centerPos, Eigen::Vector2d pointPos, double cellSize);
        double getVelCompV(Eigen::Vector2i index_Corner, Eigen::Vector2d centerPos, Eigen::Vector2d pointPos, double cellSize);

        void outputVelField();

        std::shared_ptr<SimParameters> mParams_;
        std::vector<std::vector<std::shared_ptr<Cell>>> mCells_;
        std::vector<std::shared_ptr<Cell>> mUActiveCells_;
        std::vector<std::shared_ptr<Cell>> mVActiveCells_;
        std::vector<std::shared_ptr<Cell>> mFluidCells_;
        std::vector<std::shared_ptr<MarkerParticle>> mMarkerParticles_;
        Eigen::MatrixX2d mParticleVels_;
        double mCellSize;
        double mSubTimeStep;
        //int mCellRows, mCellCols;
    };

}