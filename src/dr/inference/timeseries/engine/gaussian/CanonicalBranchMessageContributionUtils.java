package dr.inference.timeseries.engine.gaussian;

import dr.inference.timeseries.representation.CanonicalGaussianState;
import dr.inference.timeseries.representation.CanonicalGaussianUtils;

/**
 * Utilities for recovering local canonical branch contributions from canonical
 * pair states.
 */
public final class CanonicalBranchMessageContributionUtils {

    private CanonicalBranchMessageContributionUtils() {
        // no instances
    }

    public static void fillFromPairState(final CanonicalGaussianState pairState,
                                         final CanonicalBranchMessageContribution out) {
        final int pairDimension = pairState.getDimension();
        if ((pairDimension & 1) != 0) {
            throw new IllegalArgumentException("pairState dimension must be even");
        }
        final int d = pairDimension / 2;
        if (out.dLogL_dInformationX.length != d) {
            throw new IllegalArgumentException("Contribution dimension mismatch");
        }

        final double[] jointMean = new double[pairDimension];
        final double[][] jointCovariance = new double[pairDimension][pairDimension];
        CanonicalGaussianUtils.fillMomentsFromCanonical(pairState, jointMean, jointCovariance);

        for (int i = 0; i < d; ++i) {
            final double currentMeanI = jointMean[i];
            final double nextMeanI = jointMean[d + i];
            out.dLogL_dInformationX[i] = currentMeanI;
            out.dLogL_dInformationY[i] = nextMeanI;
            for (int j = 0; j < d; ++j) {
                final double currentMeanJ = jointMean[j];
                final double nextMeanJ = jointMean[d + j];
                final double currentSecondMoment = jointCovariance[i][j] + currentMeanI * currentMeanJ;
                final double crossSecondMoment = jointCovariance[d + i][j] + nextMeanI * currentMeanJ;
                final double nextSecondMoment = jointCovariance[d + i][d + j] + nextMeanI * nextMeanJ;
                out.dLogL_dPrecisionXX[i][j] = -0.5 * currentSecondMoment;
                out.dLogL_dPrecisionXY[i][j] = -0.5 * (jointCovariance[i][d + j] + currentMeanI * nextMeanJ);
                out.dLogL_dPrecisionYX[i][j] = -0.5 * crossSecondMoment;
                out.dLogL_dPrecisionYY[i][j] = -0.5 * nextSecondMoment;
            }
        }
        out.dLogL_dLogNormalizer = -1.0;
    }
}
