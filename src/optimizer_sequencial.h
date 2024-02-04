#pragma once
#include "hierarchy.h"

class Optimizer_sequencial {
public:
    Optimizer_sequencial(MultiResolutionHierarchy &mRes);

    void setOptimizeOrientations(bool value);
    void setOptimizePositions(bool value);
    void setAlignment(bool alignment) { mAlignment = alignment; }
    void setRandomization(bool randomization) { mRandomization = randomization; }
    void setExtrinsic(bool extrinsic) { mExtrinsic = extrinsic; }
    void setHierarchy(bool hierarchy) { mHierarchy = hierarchy; }
    void setMaxIterations(uint32_t it) { mMaxIterations = it; }

    bool hierarchy() const { return mHierarchy; }
    bool randomization() const { return mRandomization; }
    bool alignment() const { return mAlignment; }
    bool extrinsic() const { return mExtrinsic; }
    uint32_t maxIterations() const { return mMaxIterations; }

    void run(); 

private:
    MultiResolutionHierarchy &mRes;

    bool mRunning;
    bool mOptimizeOrientations;
    bool mOptimizePositions;
    bool mAlignment;
    bool mRandomization;
    bool mExtrinsic;
    bool mHierarchy;
    uint32_t mLevel;
    uint32_t mLevelIterations;
    uint32_t mMaxIterations;
};
