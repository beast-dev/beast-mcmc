package dr.inference.timeseries.engine.kalman.formula;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianMatrixOps;
import dr.inference.timeseries.engine.kalman.BranchSmootherStats;

/**
 * Reusable branch-local adjoint workspace for expectation-path OU gradients.
 *
 * <p>The selection, diffusion, and stationary-mean formulas all need the same
 * Gaussian transition derivatives. This helper owns that algebra once and lets
 * the callers decide how to push the adjoints back to model parameters.</p>
 */
final class ExpectationGaussianBranchGradientAdjoints {

    private final int stateDimension;

    private final double[] meanResidual;
    private final double[][] crossCov;
    private final double[][] branchMat;
    private final double[][] stepCovInv;
    private final double[][] tempDxD;
    private final double[][] tempDxD2;
    private final double[][] residualSecondMoment;
    private final double[][] dLogL_dF;
    private final double[] dLogL_df;
    private final double[][] dLogL_dV;

    ExpectationGaussianBranchGradientAdjoints(final int stateDimension) {
        if (stateDimension < 1) {
            throw new IllegalArgumentException("stateDimension must be at least 1");
        }
        this.stateDimension = stateDimension;

        this.meanResidual = new double[stateDimension];
        this.crossCov = new double[stateDimension][stateDimension];
        this.branchMat = new double[stateDimension][stateDimension];
        this.stepCovInv = new double[stateDimension][stateDimension];
        this.tempDxD = new double[stateDimension][stateDimension];
        this.tempDxD2 = new double[stateDimension][stateDimension];
        this.residualSecondMoment = new double[stateDimension][stateDimension];
        this.dLogL_dF = new double[stateDimension][stateDimension];
        this.dLogL_df = new double[stateDimension];
        this.dLogL_dV = new double[stateDimension][stateDimension];
    }

    void compute(final BranchSmootherStats curr,
                 final BranchSmootherStats next,
                 final double[][] transitionMatrix,
                 final double[] transitionOffset,
                 final double[][] stepCovariance) {
        final int d = stateDimension;

        GaussianMatrixOps.multiplyMatrixMatrixTransposedRight(
                next.smoothedCovariance, curr.smootherGain, crossCov);

        GaussianMatrixOps.multiplyMatrixVector(transitionMatrix, curr.smoothedMean, meanResidual);
        for (int i = 0; i < d; ++i) {
            meanResidual[i] = next.smoothedMean[i] - meanResidual[i] - transitionOffset[i];
        }

        GaussianMatrixOps.copyMatrix(crossCov, branchMat);
        GaussianMatrixOps.multiplyMatrixMatrix(transitionMatrix, curr.smoothedCovariance, tempDxD);
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                branchMat[i][j] -= tempDxD[i][j];
                branchMat[i][j] += meanResidual[i] * curr.smoothedMean[j];
            }
        }

        GaussianMatrixOps.copyMatrix(stepCovariance, stepCovInv);
        final GaussianMatrixOps.CholeskyFactor stepChol =
                GaussianMatrixOps.cholesky(stepCovInv);
        GaussianMatrixOps.invertPositiveDefiniteFromCholesky(stepCovInv, stepChol);

        GaussianMatrixOps.multiplyMatrixMatrix(stepCovInv, branchMat, dLogL_dF);
        GaussianMatrixOps.multiplyMatrixVector(stepCovInv, meanResidual, dLogL_df);

        GaussianMatrixOps.copyMatrix(next.smoothedCovariance, residualSecondMoment);
        GaussianMatrixOps.multiplyMatrixMatrixTransposedRight(transitionMatrix, crossCov, tempDxD);
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                residualSecondMoment[i][j] -= tempDxD[i][j];
            }
        }

        GaussianMatrixOps.multiplyMatrixMatrixTransposedRight(crossCov, transitionMatrix, tempDxD);
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                residualSecondMoment[i][j] -= tempDxD[i][j];
            }
        }

        GaussianMatrixOps.multiplyMatrixMatrix(transitionMatrix, curr.smoothedCovariance, tempDxD);
        GaussianMatrixOps.multiplyMatrixMatrixTransposedRight(tempDxD, transitionMatrix, tempDxD2);
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                residualSecondMoment[i][j] += tempDxD2[i][j];
                residualSecondMoment[i][j] += meanResidual[i] * meanResidual[j];
            }
        }

        GaussianMatrixOps.multiplyMatrixMatrix(stepCovInv, residualSecondMoment, tempDxD);
        GaussianMatrixOps.multiplyMatrixMatrix(tempDxD, stepCovInv, dLogL_dV);
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                dLogL_dV[i][j] = 0.5 * (dLogL_dV[i][j] - stepCovInv[i][j]);
            }
        }
    }

    void computeMeanAdjoint(final BranchSmootherStats curr,
                            final BranchSmootherStats next,
                            final double[][] transitionMatrix,
                            final double[] transitionOffset,
                            final double[][] stepCovariance) {
        final int d = stateDimension;

        GaussianMatrixOps.multiplyMatrixVector(transitionMatrix, curr.smoothedMean, meanResidual);
        for (int i = 0; i < d; ++i) {
            meanResidual[i] = next.smoothedMean[i] - meanResidual[i] - transitionOffset[i];
        }

        GaussianMatrixOps.copyMatrix(stepCovariance, stepCovInv);
        final GaussianMatrixOps.CholeskyFactor stepChol =
                GaussianMatrixOps.cholesky(stepCovInv);
        GaussianMatrixOps.invertPositiveDefiniteFromCholesky(stepCovInv, stepChol);
        GaussianMatrixOps.multiplyMatrixVector(stepCovInv, meanResidual, dLogL_df);
    }

    double[][] dLogL_dF() {
        return dLogL_dF;
    }

    double[] dLogL_df() {
        return dLogL_df;
    }

    double[][] dLogL_dV() {
        return dLogL_dV;
    }

    void copyDLogLDFToFlat(final double[] out) {
        copyMatrixToFlat(dLogL_dF, out);
    }

    void copyDLogLDVToFlat(final double[] out) {
        copyMatrixToFlat(dLogL_dV, out);
    }

    private void copyMatrixToFlat(final double[][] source,
                                  final double[] out) {
        for (int row = 0; row < stateDimension; ++row) {
            System.arraycopy(source[row], 0, out, row * stateDimension, stateDimension);
        }
    }
}
