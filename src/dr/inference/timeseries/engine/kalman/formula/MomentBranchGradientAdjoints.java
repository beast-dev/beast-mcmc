package dr.inference.timeseries.engine.kalman.formula;

import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.inference.timeseries.engine.kalman.MomentBranchSmootherStats;

/**
 * Reusable branch-local adjoint workspace for moment-path OU gradients.
 *
 * <p>The selection, diffusion, and stationary-mean formulas all need the same
 * Gaussian transition derivatives. This helper owns that algebra once and lets
 * the callers decide how to push the adjoints back to model parameters.</p>
 */
final class MomentBranchGradientAdjoints {

    private final int stateDimension;

    private final double[] meanResidual;
    private final double[] crossCov;
    private final double[] branchMat;
    private final double[] stepCovInv;
    private final double[] stepCovCholesky;
    private final double[] stepCovLowerInverse;
    private final double[] tempDxD;
    private final double[] tempDxD2;
    private final double[] residualSecondMoment;
    private final double[] dLogL_dF;
    private final double[] dLogL_df;
    private final double[] dLogL_dV;

    MomentBranchGradientAdjoints(final int stateDimension) {
        if (stateDimension < 1) {
            throw new IllegalArgumentException("stateDimension must be at least 1");
        }
        this.stateDimension = stateDimension;

        this.meanResidual = new double[stateDimension];
        final int matrixSize = stateDimension * stateDimension;
        this.crossCov = new double[matrixSize];
        this.branchMat = new double[matrixSize];
        this.stepCovInv = new double[matrixSize];
        this.stepCovCholesky = new double[matrixSize];
        this.stepCovLowerInverse = new double[matrixSize];
        this.tempDxD = new double[matrixSize];
        this.tempDxD2 = new double[matrixSize];
        this.residualSecondMoment = new double[matrixSize];
        this.dLogL_dF = new double[matrixSize];
        this.dLogL_df = new double[stateDimension];
        this.dLogL_dV = new double[matrixSize];
    }

    void compute(final MomentBranchSmootherStats curr,
                 final MomentBranchSmootherStats next,
                 final double[] transitionMatrix,
                 final int transitionMatrixOffset,
                 final double[] transitionOffset,
                 final int transitionOffsetOffset,
                 final double[] stepCovariance,
                 final int stepCovarianceOffset) {
        final int d = stateDimension;
        final int matrixSize = d * d;

        MatrixOps.matMulTransposedRight(
                next.smoothedCovariance, curr.smootherGain, crossCov, d);

        MatrixOps.matVec(transitionMatrix, transitionMatrixOffset, curr.smoothedMean, meanResidual, d);
        for (int i = 0; i < d; ++i) {
            meanResidual[i] = next.smoothedMean[i] - meanResidual[i] - transitionOffset[transitionOffsetOffset + i];
        }

        System.arraycopy(crossCov, 0, branchMat, 0, matrixSize);
        MatrixOps.matMul(transitionMatrix, transitionMatrixOffset, curr.smoothedCovariance, 0, tempDxD, d);
        for (int i = 0; i < d; ++i) {
            final int rowOffset = i * d;
            for (int j = 0; j < d; ++j) {
                branchMat[rowOffset + j] -= tempDxD[rowOffset + j];
                branchMat[rowOffset + j] += meanResidual[i] * curr.smoothedMean[j];
            }
        }

        System.arraycopy(stepCovariance, stepCovarianceOffset, stepCovInv, 0, matrixSize);
        if (!MatrixOps.tryCholesky(stepCovInv, stepCovCholesky, d)) {
            throw new IllegalArgumentException("Matrix is not positive definite");
        }
        MatrixOps.invertFromCholesky(stepCovCholesky, stepCovLowerInverse, stepCovInv, d);

        MatrixOps.matMul(stepCovInv, branchMat, dLogL_dF, d);
        MatrixOps.matVec(stepCovInv, meanResidual, dLogL_df, d);

        System.arraycopy(next.smoothedCovariance, 0, residualSecondMoment, 0, matrixSize);
        MatrixOps.matMulTransposedRight(transitionMatrix, transitionMatrixOffset, crossCov, 0, tempDxD, d);
        for (int i = 0; i < d; ++i) {
            final int rowOffset = i * d;
            for (int j = 0; j < d; ++j) {
                residualSecondMoment[rowOffset + j] -= tempDxD[rowOffset + j];
            }
        }

        MatrixOps.matMulTransposedRight(crossCov, 0, transitionMatrix, transitionMatrixOffset, tempDxD, d);
        for (int i = 0; i < d; ++i) {
            final int rowOffset = i * d;
            for (int j = 0; j < d; ++j) {
                residualSecondMoment[rowOffset + j] -= tempDxD[rowOffset + j];
            }
        }

        MatrixOps.matMul(transitionMatrix, transitionMatrixOffset, curr.smoothedCovariance, 0, tempDxD, d);
        MatrixOps.matMulTransposedRight(tempDxD, 0, transitionMatrix, transitionMatrixOffset, tempDxD2, d);
        for (int i = 0; i < d; ++i) {
            final int rowOffset = i * d;
            for (int j = 0; j < d; ++j) {
                residualSecondMoment[rowOffset + j] += tempDxD2[rowOffset + j];
                residualSecondMoment[rowOffset + j] += meanResidual[i] * meanResidual[j];
            }
        }

        MatrixOps.matMul(stepCovInv, residualSecondMoment, tempDxD, d);
        MatrixOps.matMul(tempDxD, stepCovInv, dLogL_dV, d);
        for (int i = 0; i < d; ++i) {
            final int rowOffset = i * d;
            for (int j = 0; j < d; ++j) {
                dLogL_dV[rowOffset + j] = 0.5 * (dLogL_dV[rowOffset + j] - stepCovInv[rowOffset + j]);
            }
        }
    }

    void computeMeanAdjoint(final MomentBranchSmootherStats curr,
                            final MomentBranchSmootherStats next,
                            final double[] transitionMatrix,
                            final int transitionMatrixOffset,
                            final double[] transitionOffset,
                            final int transitionOffsetOffset,
                            final double[] stepCovariance,
                            final int stepCovarianceOffset) {
        final int d = stateDimension;
        final int matrixSize = d * d;

        MatrixOps.matVec(transitionMatrix, transitionMatrixOffset, curr.smoothedMean, meanResidual, d);
        for (int i = 0; i < d; ++i) {
            meanResidual[i] = next.smoothedMean[i] - meanResidual[i] - transitionOffset[transitionOffsetOffset + i];
        }

        System.arraycopy(stepCovariance, stepCovarianceOffset, stepCovInv, 0, matrixSize);
        if (!MatrixOps.tryCholesky(stepCovInv, stepCovCholesky, d)) {
            throw new IllegalArgumentException("Matrix is not positive definite");
        }
        MatrixOps.invertFromCholesky(stepCovCholesky, stepCovLowerInverse, stepCovInv, d);
        MatrixOps.matVec(stepCovInv, meanResidual, dLogL_df, d);
    }

    double[] dLogL_dF() {
        return dLogL_dF;
    }

    double[] dLogL_df() {
        return dLogL_df;
    }

    double[] dLogL_dV() {
        return dLogL_dV;
    }

    void copyDLogLDFToFlat(final double[] out) {
        System.arraycopy(dLogL_dF, 0, out, 0, dLogL_dF.length);
    }

    void copyDLogLDVToFlat(final double[] out) {
        System.arraycopy(dLogL_dV, 0, out, 0, dLogL_dV.length);
    }
}
