#include "optimizer_sequencial.h"

Optimizer_sequencial::Optimizer_sequencial(MultiResolutionHierarchy &mRes)
    : mRes(mRes), mRunning(false), mOptimizeOrientations(false),
      mOptimizePositions(false), mAlignment(true), mRandomization(true), mExtrinsic(true),
      mHierarchy(true), mLevel(0), mLevelIterations(0), mMaxIterations(20) {
    ;
}

void Optimizer_sequencial::setOptimizeOrientations(bool value) {
    std::lock_guard<ordered_lock> lock(mRes.mutex());
	if (mRes.D3) mMaxIterations = 200;else mMaxIterations = 20;
    mOptimizeOrientations = value;
    mLevel = mRes.levels() - 2;
    mLevelIterations = 0;
}

void Optimizer_sequencial::setOptimizePositions(bool value) {
    std::lock_guard<ordered_lock> lock(mRes.mutex());
	if (mRes.D3) mMaxIterations = 200; else mMaxIterations = 50;
	mOptimizePositions = value;
    mLevel = mRes.levels() - 2;
    mLevelIterations = 0;
}

void Optimizer_sequencial::run() {
    mRunning = true;

    while (true) {
        if (!mHierarchy)
            mLevel = 0;

        if (mOptimizeOrientations) {
            if (mRes.D3)
                mRes.smoothOrientationsTet(mLevel, mAlignment, mRandomization);
            else
                mRes.smoothOrientationsTri(mLevel, mAlignment, mRandomization, mExtrinsic);
        }

        if (mOptimizePositions) {
            if (mRes.D3)
                mRes.smoothPositionsTet(mLevel, mAlignment, mRandomization);
            else {
                mRes.smoothPositionsTri(mLevel, mAlignment, mRandomization, mExtrinsic);
            }
        }

        mLevelIterations++;

		if (mLevelIterations >= mMaxIterations) {
			mLevelIterations = 0;
			if (mLevel == 0) {
				mOptimizeOrientations = false;
				mOptimizePositions = false;
				break;
			}
			if (mHierarchy) {
				mLevel--;
				if (mOptimizeOrientations)
					mRes.prolongOrientations(mLevel);
				if (mOptimizePositions)
					mRes.prolongPositions(mLevel);
			}
		}
    }
}