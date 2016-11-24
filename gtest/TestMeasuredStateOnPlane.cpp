#include <gtest/gtest.h>

#include <SharedPlanePtr.h>
#include <RKTrackRep.h>
#include <MeasuredStateOnPlane.h>

#include <TMatrixD.h>
#include <TVector3.h>

namespace genfit {

    class MeasuredStateOnPlaneTests : public ::testing::Test {
    protected:
        static constexpr double numericalTolerance = 1e-10;
    };

    TEST_F(MeasuredStateOnPlaneTests, Test) {
        RKTrackRep myRKTrackRep;
        genfit::SharedPlanePtr mySharedPlanePtr(new genfit::DetPlane(TVector3(0, 0, 0), TVector3(0, 0, 1), nullptr));

        TVectorD state1(5);
        state1(0) = 1;
        state1(1) = 0.1;
        state1(2) = 0.1;
        state1(3) = 0.01;
        state1(4) = 0.01;

        TVectorD state2(5);
        state2(0) = 1;
        state2(1) = -0.1;
        state2(2) = -0.1;
        state2(3) = -0.01;
        state2(4) = -0.01;

        TMatrixDSym covariance1(5);
        TMatrixDSym covariance2(5);
        for (unsigned int row=0; row<5; ++row) {
            for (unsigned int col=0; col<5; ++col) {
                if (row == col) {
                    covariance1(row, col) = 0.1;
                    covariance2(row, col) = 0.1;
                } else {
                    covariance1(row, col) = 0.01;
                    covariance2(row, col) = 0.01;
                }
            }
        }

        genfit::MeasuredStateOnPlane myMeasuredStateOnPlane1(
                state1, covariance1, mySharedPlanePtr, &myRKTrackRep);
        genfit::MeasuredStateOnPlane myMeasuredStateOnPlane2(
                state2, covariance2, mySharedPlanePtr, &myRKTrackRep);


        genfit::MeasuredStateOnPlane myAveragedStateOnPlane(calcAverageState(myMeasuredStateOnPlane1, myMeasuredStateOnPlane2));
        const auto myAveragedState = myAveragedStateOnPlane.getState();
        const auto myAveragedCov = myAveragedStateOnPlane.getCov();

        EXPECT_NEAR(1, myAveragedState(0), numericalTolerance);
        EXPECT_NEAR(0, myAveragedState(1), numericalTolerance);
        EXPECT_NEAR(0, myAveragedState(2), numericalTolerance);
        EXPECT_NEAR(0, myAveragedState(3), numericalTolerance);
        EXPECT_NEAR(0, myAveragedState(4), numericalTolerance);

        // diagonal elements
        EXPECT_NEAR(0.05, myAveragedCov(0, 0), numericalTolerance);
        EXPECT_NEAR(0.05, myAveragedCov(1, 1), numericalTolerance);
        EXPECT_NEAR(0.05, myAveragedCov(2, 2), numericalTolerance);
        EXPECT_NEAR(0.05, myAveragedCov(3, 3), numericalTolerance);
        EXPECT_NEAR(0.05, myAveragedCov(4, 4), numericalTolerance);

        // off-diagonal elements
        EXPECT_NEAR(0.005, myAveragedCov(0, 1), numericalTolerance);
        EXPECT_NEAR(0.005, myAveragedCov(0, 2), numericalTolerance);
        EXPECT_NEAR(0.005, myAveragedCov(0, 3), numericalTolerance);
        EXPECT_NEAR(0.005, myAveragedCov(0, 4), numericalTolerance);
        EXPECT_NEAR(0.005, myAveragedCov(1, 2), numericalTolerance);
        EXPECT_NEAR(0.005, myAveragedCov(1, 3), numericalTolerance);
        EXPECT_NEAR(0.005, myAveragedCov(1, 4), numericalTolerance);
        EXPECT_NEAR(0.005, myAveragedCov(2, 3), numericalTolerance);
        EXPECT_NEAR(0.005, myAveragedCov(2, 4), numericalTolerance);
        EXPECT_NEAR(0.005, myAveragedCov(3, 4), numericalTolerance);

    }

}
