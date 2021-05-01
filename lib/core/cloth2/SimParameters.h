#ifndef PSIM_CORE_CLOTH2_SIMPARAMETERS_H
#define PSIM_CORE_CLOTH2_SIMPARAMETERS_H

namespace cloth2 {

struct SimParameters {
    float dt = 1e-3;
    int constraintIters = 5;

    bool gravityEnabled = true;
    float gravityG = -9.8;

    bool pinEnabled = true;
    float pinWeight = 1.0;

    bool stretchEnabled = true;
    float stretchWeight = 0.5;

    bool bendingEnabled = true;
    float bendingWeight = 0.5;

    bool dragEnabled = true;
    float dragWeight = 0.5;

    bool collisionsEnabled = true;

    bool shouldDeployBevo = false;
    float launchSpeed = 100.0;
    float launchSpinningSpeed = 25.0;
    
    bool enableIntermediateValueRecording = false;
    bool enableCVBF = true;
    bool enableCFBV = true;
};

}

#endif
