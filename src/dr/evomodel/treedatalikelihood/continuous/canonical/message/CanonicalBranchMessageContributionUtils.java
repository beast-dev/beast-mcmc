package dr.evomodel.treedatalikelihood.continuous.canonical.message;

import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;

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
        final double[] precisionXX;
        final double[] precisionXY;
        final double[] precisionYX;
        final double[] precisionYY;
        final double[] inverseXX;
        final double[] inverseXXTimesXY;
        final double[] schurYY;
        final double[] covarianceYY;
        final double[] covarianceXX;
        final double[] covarianceXY;
        final double[] matrixScratch;
        final double[] choleskyScratch;
        final double[] lowerInverseScratch;
        final double[] symmetricScratch;
        private final int stateDimension;

        public Workspace(final int stateDimension) {
            this.stateDimension = stateDimension;
            final int pairDimension = 2 * stateDimension;
            final int d2 = stateDimension * stateDimension;
            this.jointMean = new double[pairDimension];
            this.precisionXX = new double[d2];
            this.precisionXY = new double[d2];
            this.precisionYX = new double[d2];
            this.precisionYY = new double[d2];
            this.inverseXX = new double[d2];
            this.inverseXXTimesXY = new double[d2];
            this.schurYY = new double[d2];
            this.covarianceYY = new double[d2];
            this.covarianceXX = new double[d2];
            this.covarianceXY = new double[d2];
            this.matrixScratch = new double[d2];
            this.choleskyScratch = new double[d2];
            this.lowerInverseScratch = new double[d2];
            this.symmetricScratch = new double[d2];
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

        fillPairMomentsFromCanonicalBlocks(pairState, workspace, d);

        final double[] jointMean = workspace.jointMean;
        final double[] covarianceXX = workspace.covarianceXX;
        final double[] covarianceXY = workspace.covarianceXY;
        final double[] covarianceYY = workspace.covarianceYY;
        for (int i = 0; i < d; ++i) {
            final double currentMeanI = jointMean[i];
            final double nextMeanI = jointMean[d + i];
            out.dLogL_dInformationX[i] = currentMeanI;
            out.dLogL_dInformationY[i] = nextMeanI;
            for (int j = 0; j < d; ++j) {
                final double currentMeanJ = jointMean[j];
                final double nextMeanJ = jointMean[d + j];
                final int ij = i * d + j;
                final double currentSecondMoment =
                        covarianceXX[ij] + currentMeanI * currentMeanJ;
                final double crossSecondMoment =
                        covarianceXY[j * d + i] + nextMeanI * currentMeanJ;
                final double nextSecondMoment =
                        covarianceYY[ij] + nextMeanI * nextMeanJ;
                out.dLogL_dPrecisionXX[ij] = -0.5 * currentSecondMoment;
                out.dLogL_dPrecisionXY[ij] =
                        -0.5 * (covarianceXY[ij] + currentMeanI * nextMeanJ);
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

    private static void fillPairMomentsFromCanonicalBlocks(final CanonicalGaussianState pairState,
                                                           final Workspace workspace,
                                                           final int d) {
        final int pairDimension = 2 * d;
        final double[] precision = pairState.precision;

        for (int i = 0; i < d; ++i) {
            final int upperRowOffset = i * pairDimension;
            final int lowerRowOffset = (d + i) * pairDimension;
            final int blockRowOffset = i * d;
            for (int j = 0; j < d; ++j) {
                workspace.precisionXX[blockRowOffset + j] = precision[upperRowOffset + j];
                workspace.precisionXY[blockRowOffset + j] = precision[upperRowOffset + d + j];
                workspace.precisionYX[blockRowOffset + j] = precision[lowerRowOffset + j];
                workspace.precisionYY[blockRowOffset + j] = precision[lowerRowOffset + d + j];
            }
        }

        invertPositiveDefiniteSymmetric(
                workspace.precisionXX,
                workspace.inverseXX,
                d,
                workspace);
        MatrixOps.matMul(
                workspace.inverseXX,
                workspace.precisionXY,
                workspace.inverseXXTimesXY,
                d);
        MatrixOps.matMul(
                workspace.precisionYX,
                workspace.inverseXXTimesXY,
                workspace.matrixScratch,
                d);
        for (int k = 0; k < d * d; ++k) {
            workspace.schurYY[k] = workspace.precisionYY[k] - workspace.matrixScratch[k];
        }
        MatrixOps.symmetrize(workspace.schurYY, d);

        invertPositiveDefiniteSymmetric(
                workspace.schurYY,
                workspace.covarianceYY,
                d,
                workspace);

        MatrixOps.matMul(
                workspace.inverseXXTimesXY,
                workspace.covarianceYY,
                workspace.matrixScratch,
                d);
        MatrixOps.copyMatrix(workspace.matrixScratch, workspace.covarianceXY, d);
        MatrixOps.scaleInPlace(workspace.covarianceXY, -1.0, d * d);

        MatrixOps.matMulTransposedRight(
                workspace.matrixScratch,
                workspace.inverseXXTimesXY,
                workspace.covarianceXX,
                d);
        MatrixOps.addInPlace(workspace.covarianceXX, workspace.inverseXX, d * d);
        MatrixOps.symmetrize(workspace.covarianceXX, d);

        final double[] h = pairState.information;
        final double[] mean = workspace.jointMean;
        for (int i = 0; i < d; ++i) {
            double meanX = 0.0;
            double meanY = 0.0;
            final int rowOffset = i * d;
            for (int j = 0; j < d; ++j) {
                meanX += workspace.covarianceXX[rowOffset + j] * h[j]
                        + workspace.covarianceXY[rowOffset + j] * h[d + j];
                meanY += workspace.covarianceXY[j * d + i] * h[j]
                        + workspace.covarianceYY[rowOffset + j] * h[d + j];
            }
            mean[i] = meanX;
            mean[d + i] = meanY;
        }
    }

    private static void invertPositiveDefiniteSymmetric(final double[] matrix,
                                                        final double[] inverseOut,
                                                        final int d,
                                                        final Workspace workspace) {
        for (int i = 0; i < d; ++i) {
            final int rowOffset = i * d;
            for (int j = 0; j < d; ++j) {
                workspace.symmetricScratch[rowOffset + j] =
                        0.5 * (matrix[rowOffset + j] + matrix[j * d + i]);
            }
        }
        if (!MatrixOps.tryCholesky(workspace.symmetricScratch, workspace.choleskyScratch, d)) {
            throw new IllegalArgumentException("Pair precision block is not positive definite.");
        }
        MatrixOps.invertFromCholesky(
                workspace.choleskyScratch,
                workspace.lowerInverseScratch,
                inverseOut,
                d);
    }
}
