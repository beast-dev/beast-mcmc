package dr.evomodel.treedatalikelihood.continuous.canonical.message;

import dr.evomodel.treedatalikelihood.continuous.canonical.math.GaussianFormConverter;

/**
 * Utilities for recovering local canonical branch contributions from canonical
 * pair states.
 *
 * <p>Callers are expected to retain and reuse a {@link Workspace}; the hot-path
 * method in this helper does not allocate.
 */
public final class CanonicalBranchMessageContributionUtils {

    public static final class Workspace {
        final double[] jointMean;
        final double[] jointCovarianceFlat;
        final GaussianFormConverter.Workspace converterWorkspace;
        private final int stateDimension;

        public Workspace(final int stateDimension) {
            this.stateDimension = stateDimension;
            final int pairDimension = 2 * stateDimension;
            this.jointMean = new double[pairDimension];
            this.jointCovarianceFlat = new double[pairDimension * pairDimension];
            this.converterWorkspace = new GaussianFormConverter.Workspace();
            this.converterWorkspace.ensureDim(pairDimension);
        }

        public int getStateDimension() {
            return stateDimension;
        }
    }

    private CanonicalBranchMessageContributionUtils() {
        // no instances
    }

    public static void fillFromPairState(final CanonicalGaussianState pairState,
                                         final Workspace workspace,
                                         final CanonicalBranchMessageContribution out) {
        final int pairDimension = pairState.getDimension();
        if ((pairDimension & 1) != 0) {
            throw new IllegalArgumentException("pairState dimension must be even");
        }
        final int d = pairDimension / 2;
        ensureWorkspaceDimension(workspace, d);
        if (out.getDimension() != d) {
            throw new IllegalArgumentException("Contribution dimension mismatch");
        }

        final double[] jointMean = workspace.jointMean;
        final double[] jointCovariance = workspace.jointCovarianceFlat;
        GaussianFormConverter.fillMomentsFromState(
                pairState,
                jointMean,
                workspace.jointCovarianceFlat,
                pairDimension,
                workspace.converterWorkspace);

        for (int i = 0; i < d; ++i) {
            final double currentMeanI = jointMean[i];
            final double nextMeanI = jointMean[d + i];
            out.dLogL_dInformationX[i] = currentMeanI;
            out.dLogL_dInformationY[i] = nextMeanI;
            for (int j = 0; j < d; ++j) {
                final double currentMeanJ = jointMean[j];
                final double nextMeanJ = jointMean[d + j];
                final double currentSecondMoment =
                        jointCovariance[i * pairDimension + j] + currentMeanI * currentMeanJ;
                final double crossSecondMoment =
                        jointCovariance[(d + i) * pairDimension + j] + nextMeanI * currentMeanJ;
                final double nextSecondMoment =
                        jointCovariance[(d + i) * pairDimension + d + j] + nextMeanI * nextMeanJ;
                final int ij = i * d + j;
                out.dLogL_dPrecisionXX[ij] = -0.5 * currentSecondMoment;
                out.dLogL_dPrecisionXY[ij] =
                        -0.5 * (jointCovariance[i * pairDimension + d + j] + currentMeanI * nextMeanJ);
                out.dLogL_dPrecisionYX[ij] = -0.5 * crossSecondMoment;
                out.dLogL_dPrecisionYY[ij] = -0.5 * nextSecondMoment;
            }
        }
        out.dLogL_dLogNormalizer = -1.0;
    }

    private static void ensureWorkspaceDimension(final Workspace workspace, final int stateDimension) {
        if (workspace.getStateDimension() != stateDimension) {
            throw new IllegalArgumentException("Workspace dimension mismatch");
        }
    }
}
