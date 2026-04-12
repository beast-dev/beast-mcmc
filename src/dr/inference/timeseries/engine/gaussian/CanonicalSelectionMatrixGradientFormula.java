package dr.inference.timeseries.engine.gaussian;

import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.OrthogonalBlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.OrthogonalMatrixProvider;
import dr.inference.model.Parameter;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.gaussian.OUProcessModel;
import dr.inference.timeseries.gaussian.OrthogonalBlockDiagonalSelectionMatrixParameterization;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;

/**
 * Selection-matrix gradient formula driven by branch-local canonical posterior states.
 *
 * <p>This implementation is the first step toward native canonical branch adjoints:
 * for each branch it constructs the joint posterior over {@code (x_t, x_{t+1})} in
 * canonical form, then derives the branch moments needed to build the usual
 * transition-space sensitivities {@code dL/dF}, {@code dL/df}, and {@code dL/dV}.
 */
public final class CanonicalSelectionMatrixGradientFormula implements CanonicalGradientFormula {

    private final Parameter selectionMatrixParameter;
    private final int stateDimension;
    private final OUProcessModel processModel;

    private final CanonicalBranchPosterior branchPosterior;
    private final double[][] stepCovInv;
    private final double[][] tempDxD;
    private final CanonicalBranchMessageContribution localContribution;
    private final CanonicalLocalTransitionAdjoints localAdjoints;
    private final double[][] branchTransitionMatrix;
    private final double[] branchTransitionOffset;
    private final double[][] tempDxD2;
    private final double[] dLogL_df;
    private final double[][] dLogL_dF;
    private final double[][] dLogL_dV;
    private final double[] nativeCompressedGradientScratch;
    private final double[][] nativeRotationGradientScratch;
    private final double[] stationaryMeanScratch;

    public CanonicalSelectionMatrixGradientFormula(final Parameter selectionMatrixParameter,
                                                   final int stateDimension) {
        this(null, selectionMatrixParameter, stateDimension);
    }

    public CanonicalSelectionMatrixGradientFormula(final OUProcessModel processModel,
                                                   final Parameter selectionMatrixParameter,
                                                   final int stateDimension) {
        if (selectionMatrixParameter == null) {
            throw new IllegalArgumentException("selectionMatrixParameter must not be null");
        }
        if (stateDimension < 1) {
            throw new IllegalArgumentException("stateDimension must be at least 1");
        }
        this.selectionMatrixParameter = selectionMatrixParameter;
        this.stateDimension = stateDimension;
        this.processModel = processModel;

        this.branchPosterior = new CanonicalBranchPosterior(stateDimension);
        this.stepCovInv = new double[stateDimension][stateDimension];
        this.tempDxD = new double[stateDimension][stateDimension];
        this.localContribution = new CanonicalBranchMessageContribution(stateDimension);
        this.localAdjoints = new CanonicalLocalTransitionAdjoints(stateDimension);
        this.branchTransitionMatrix = new double[stateDimension][stateDimension];
        this.branchTransitionOffset = new double[stateDimension];
        this.tempDxD2 = new double[stateDimension][stateDimension];
        this.dLogL_df = new double[stateDimension];
        this.dLogL_dF = new double[stateDimension][stateDimension];
        this.dLogL_dV = new double[stateDimension][stateDimension];
        this.nativeCompressedGradientScratch = selectionMatrixParameter instanceof AbstractBlockDiagonalTwoByTwoMatrixParameter
                ? new double[((AbstractBlockDiagonalTwoByTwoMatrixParameter) selectionMatrixParameter).getCompressedDDimension()]
                : null;
        this.nativeRotationGradientScratch = selectionMatrixParameter instanceof AbstractBlockDiagonalTwoByTwoMatrixParameter
                ? new double[stateDimension][stateDimension]
                : null;
        this.stationaryMeanScratch = processModel != null ? new double[stateDimension] : null;
    }

    @Override
    public boolean supportsParameter(final Parameter parameter) {
        if (parameter == selectionMatrixParameter) {
            return true;
        }
        if (!(selectionMatrixParameter instanceof AbstractBlockDiagonalTwoByTwoMatrixParameter)) {
            return false;
        }

        final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter =
                (AbstractBlockDiagonalTwoByTwoMatrixParameter) selectionMatrixParameter;
        if (parameter == blockParameter.getParameter()) {
            return true;
        }
        if (parameter == blockParameter.getRotationMatrixParameter()) {
            return true;
        }
        if (blockParameter.getRotationMatrixParameter() instanceof OrthogonalMatrixProvider
                && parameter == ((OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter()).getOrthogonalParameter()) {
            return true;
        }
        if (parameter == blockParameter.getScalarBlockParameter()
                && blockParameter.getScalarBlockParameter().getDimension() > 0) {
            return true;
        }
        for (int i = 0; i < blockParameter.getTwoByTwoParameterFamilyCount(); ++i) {
            if (parameter == blockParameter.getTwoByTwoBlockParameter(i)) {
                return true;
            }
        }
        return false;
    }

    @Override
    public double[] computeGradient(final Parameter parameter,
                                    final CanonicalForwardTrajectory trajectory,
                                    final GaussianTransitionRepresentation repr,
                                    final TimeGrid timeGrid) {
        if (shouldUseOrthogonalNativePath(parameter)) {
            return computeOrthogonalNativeGradient(parameter, trajectory, timeGrid);
        }

        final int d = stateDimension;
        final int T = trajectory.timeCount;
        final double[] gradientAccumulator = new double[d * d];

        for (int t = 0; t < T - 1; ++t) {
            branchPosterior.fillFromCanonicalPairState(trajectory.branchPairStates[t]);
            branchPosterior.fillLocalMessageContribution(localContribution);
            fillTransitionMoments(trajectory.transitions[t]);
            fillSelectionTransitionAdjoints(trajectory.transitions[t]);
            repr.accumulateSelectionGradient(t, t + 1, timeGrid, dLogL_dF, dLogL_df, gradientAccumulator);

            fillCovarianceAdjoint(trajectory.transitions[t]);
            repr.accumulateSelectionGradientFromCovariance(t, t + 1, timeGrid, dLogL_dV, gradientAccumulator);
        }

        return gradientAccumulator;
    }

    private boolean shouldUseOrthogonalNativePath(final Parameter requestedParameter) {
        if (requestedParameter == selectionMatrixParameter) {
            return false;
        }
        if (!(selectionMatrixParameter instanceof OrthogonalBlockDiagonalPolarStableMatrixParameter)) {
            return false;
        }
        return processModel != null
                && processModel.getCovarianceGradientMethod() == OUProcessModel.CovarianceGradientMethod.STATIONARY_LYAPUNOV
                && processModel.getSelectionMatrixParameterization()
                instanceof OrthogonalBlockDiagonalSelectionMatrixParameterization;
    }

    private double[] computeOrthogonalNativeGradient(final Parameter requestedParameter,
                                                     final CanonicalForwardTrajectory trajectory,
                                                     final TimeGrid timeGrid) {
        final OrthogonalBlockDiagonalSelectionMatrixParameterization orthogonalParameterization =
                (OrthogonalBlockDiagonalSelectionMatrixParameterization) processModel.getSelectionMatrixParameterization();
        final OrthogonalBlockDiagonalPolarStableMatrixParameter blockParameter =
                (OrthogonalBlockDiagonalPolarStableMatrixParameter) selectionMatrixParameter;
        final int T = trajectory.timeCount;
        final double[] compressedDGradient = nativeCompressedGradientScratch;
        final double[][] rotationGradient = nativeRotationGradientScratch;
        final double[] stationaryMean = stationaryMeanScratch;
        zero(compressedDGradient);
        zero(rotationGradient);
        processModel.getInitialMean(stationaryMean);

        for (int t = 0; t < T - 1; ++t) {
            final double dt = timeGrid.getDelta(t, t + 1);
            CanonicalBranchMessageContributionUtils.fillFromPairState(
                    trajectory.branchPairStates[t],
                    localContribution);
            orthogonalParameterization.accumulateNativeGradientFromCanonicalContribution(
                    processModel.getDiffusionMatrix(),
                    stationaryMean,
                    dt,
                    localContribution,
                    localAdjoints,
                    compressedDGradient, rotationGradient);
        }

        final double[] nativeGradient = new double[blockParameter.getBlockDiagonalNParameters()];
        blockParameter.chainGradient(compressedDGradient, nativeGradient);
        return assembleBlockGradientResult(requestedParameter, blockParameter, nativeGradient, rotationGradient);
    }

    private void fillTransitionMoments(final dr.inference.timeseries.representation.CanonicalGaussianTransition transition) {
        KalmanLikelihoodEngine.copyMatrix(transition.precisionYY, stepCovInv);
        final KalmanLikelihoodEngine.CholeskyFactor stepChol =
                KalmanLikelihoodEngine.cholesky(stepCovInv);
        KalmanLikelihoodEngine.invertPositiveDefiniteFromCholesky(stepCovInv, stepChol);

        KalmanLikelihoodEngine.multiplyMatrixMatrix(stepCovInv, transition.precisionYX, branchTransitionMatrix);
        for (int i = 0; i < stateDimension; ++i) {
            for (int j = 0; j < stateDimension; ++j) {
                branchTransitionMatrix[i][j] = -branchTransitionMatrix[i][j];
            }
        }
        KalmanLikelihoodEngine.multiplyMatrixVector(stepCovInv, transition.informationY, branchTransitionOffset);
    }

    private void fillSelectionTransitionAdjoints(final dr.inference.timeseries.representation.CanonicalGaussianTransition transition) {
        final int d = stateDimension;
        final double[][] gXx = localContribution.dLogL_dPrecisionXX;
        final double[][] gXy = localContribution.dLogL_dPrecisionXY;
        final double[][] gYx = localContribution.dLogL_dPrecisionYX;
        final double[] gX = localContribution.dLogL_dInformationX;
        final double[] gY = localContribution.dLogL_dInformationY;

        for (int i = 0; i < d; ++i) {
            double sum = 0.0;
            for (int k = 0; k < d; ++k) {
                sum += transition.precisionYX[i][k] * gX[k];
                sum += transition.precisionYY[i][k] * gY[k];
            }
            dLogL_df[i] = sum + localContribution.dLogL_dLogNormalizer * transition.informationY[i];
        }

        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double sum = -transition.informationY[i] * gX[j];
                for (int k = 0; k < d; ++k) {
                    sum -= transition.precisionYX[i][k] * gXx[k][j];
                    sum -= transition.precisionYX[i][k] * gXx[j][k];
                    sum -= transition.precisionYY[i][k] * gXy[j][k];
                    sum -= transition.precisionYY[i][k] * gYx[k][j];
                }
                dLogL_dF[i][j] = sum;
            }
        }
    }

    private void fillCovarianceAdjoint(final dr.inference.timeseries.representation.CanonicalGaussianTransition transition) {
        final int d = stateDimension;
        final double[][] gXx = localContribution.dLogL_dPrecisionXX;
        final double[][] gXy = localContribution.dLogL_dPrecisionXY;
        final double[][] gYx = localContribution.dLogL_dPrecisionYX;
        final double[][] gYy = localContribution.dLogL_dPrecisionYY;
        final double[] gX = localContribution.dLogL_dInformationX;
        final double[] gY = localContribution.dLogL_dInformationY;
        final double g0 = localContribution.dLogL_dLogNormalizer;

        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double sum = gYy[i][j]
                        + gY[i] * branchTransitionOffset[j]
                        + 0.5 * g0 * (branchTransitionOffset[i] * branchTransitionOffset[j] - stepCovInv[i][j]);
                for (int a = 0; a < d; ++a) {
                    sum -= branchTransitionMatrix[i][a] * gXy[a][j];
                    sum -= gYx[i][a] * branchTransitionMatrix[j][a];
                    sum -= branchTransitionMatrix[i][a] * gX[a] * branchTransitionOffset[j];
                    for (int b = 0; b < d; ++b) {
                        sum += branchTransitionMatrix[i][a] * gXx[a][b] * branchTransitionMatrix[j][b];
                    }
                }
                tempDxD[i][j] = sum;
            }
        }

        KalmanLikelihoodEngine.multiplyMatrixMatrix(transition.precisionYY, tempDxD, tempDxD2);
        KalmanLikelihoodEngine.multiplyMatrixMatrix(tempDxD2, transition.precisionYY, dLogL_dV);
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                dLogL_dV[i][j] = -dLogL_dV[i][j];
            }
        }
        KalmanLikelihoodEngine.symmetrize(dLogL_dV);
    }

    private double[] assembleBlockGradientResult(final Parameter requestedParameter,
                                                 final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                                 final double[] nativeGradient,
                                                 final double[][] gradientR) {
        final int d = stateDimension;
        if (requestedParameter == blockParameter.getParameter()) {
            if (blockParameter.getRotationMatrixParameter() instanceof OrthogonalMatrixProvider) {
                final double[] angleGradient =
                        ((OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter()).pullBackGradient(gradientR);
                final double[] out = new double[nativeGradient.length + angleGradient.length];
                System.arraycopy(nativeGradient, 0, out, 0, nativeGradient.length);
                System.arraycopy(angleGradient, 0, out, nativeGradient.length, angleGradient.length);
                return out;
            }
            final double[] out = new double[nativeGradient.length + d * d];
            System.arraycopy(nativeGradient, 0, out, 0, nativeGradient.length);
            flattenColumnMajor(gradientR, out, nativeGradient.length);
            return out;
        }
        if (requestedParameter == blockParameter.getRotationMatrixParameter()) {
            return flattenColumnMajor(gradientR);
        }
        if (blockParameter.getRotationMatrixParameter() instanceof OrthogonalMatrixProvider
                && requestedParameter == ((OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter()).getOrthogonalParameter()) {
            return ((OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter())
                    .pullBackGradient(gradientR);
        }
        if (requestedParameter == blockParameter.getScalarBlockParameter()) {
            return new double[]{nativeGradient[0]};
        }

        final int blockBase = blockParameter.hasLeadingOneByOneBlock() ? 1 : 0;
        final int blockWidth = blockParameter.getNum2x2Blocks();
        for (int family = 0; family < blockParameter.getTwoByTwoParameterFamilyCount(); ++family) {
            if (requestedParameter == blockParameter.getTwoByTwoBlockParameter(family)) {
                final double[] out = new double[blockWidth];
                System.arraycopy(nativeGradient, blockBase + family * blockWidth, out, 0, blockWidth);
                return out;
            }
        }
        throw new IllegalArgumentException("Unsupported block parameter: " + requestedParameter.getId());
    }

    private static double[] flattenColumnMajor(final double[][] matrix) {
        final double[] out = new double[matrix.length * matrix.length];
        flattenColumnMajor(matrix, out, 0);
        return out;
    }

    private static void zero(final double[] vector) {
        for (int i = 0; i < vector.length; ++i) {
            vector[i] = 0.0;
        }
    }

    private static void zero(final double[][] matrix) {
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = 0; j < matrix[i].length; ++j) {
                matrix[i][j] = 0.0;
            }
        }
    }

    private static void flattenColumnMajor(final double[][] matrix,
                                           final double[] out,
                                           final int offset) {
        final int d = matrix.length;
        int index = offset;
        for (int col = 0; col < d; ++col) {
            for (int row = 0; row < d; ++row) {
                out[index++] = matrix[row][col];
            }
        }
    }
}
