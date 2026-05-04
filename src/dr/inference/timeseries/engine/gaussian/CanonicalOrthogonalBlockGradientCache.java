package dr.inference.timeseries.engine.gaussian;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.canonical.CanonicalBranchWorkspace;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedBranchHandle;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockCanonicalParameterization;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.GaussianFormConverter;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianMatrixOps;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.OrthogonalBlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.OrthogonalMatrixProvider;
import dr.inference.model.Parameter;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.gaussian.DiffusionMatrixParameterization;

/**
 * Timeseries-only joint accumulator for orthogonal-block canonical OU gradients.
 */
public final class CanonicalOrthogonalBlockGradientCache {

    private final OUProcessModel processModel;
    private final OrthogonalBlockDiagonalPolarStableMatrixParameter blockParameter;
    private final OrthogonalBlockCanonicalParameterization orthogonalParameterization;
    private final DiffusionMatrixParameterization diffusionParameterization;
    private final Parameter stationaryMeanParameter;
    private final MatrixParameter initialCovarianceParameter;
    private final int stateDimension;

    private final double[] stationaryMean;
    private final double[] compressedDGradient;
    private final double[] rotationGradientFlat;
    private final double[][] rotationGradient2D;
    private final double[] nativeGradient;
    private final double[] diffusionGradient;
    private final double[] meanGradient;
    private final double[] smoothedInitialMean;
    private final double[] stateDiff;
    private final double[][] initialCovInv;
    private final double[][] smoothedInitialCovariance;
    private final double[] initialGradient;
    private final CanonicalPreparedBranchHandle preparedBasis;
    private final CanonicalBranchWorkspace branchWorkspace;
    private boolean known;

    public CanonicalOrthogonalBlockGradientCache(final OUProcessModel processModel,
                                                 final OrthogonalBlockDiagonalPolarStableMatrixParameter blockParameter,
                                                 final DiffusionMatrixParameterization diffusionParameterization,
                                                 final Parameter stationaryMeanParameter,
                                                 final MatrixParameter initialCovarianceParameter,
                                                 final int stateDimension) {
        this.processModel = processModel;
        this.blockParameter = blockParameter;
        this.orthogonalParameterization =
                (OrthogonalBlockCanonicalParameterization) processModel.getSelectionMatrixParameterization();
        this.diffusionParameterization = diffusionParameterization;
        this.stationaryMeanParameter = stationaryMeanParameter;
        this.initialCovarianceParameter = initialCovarianceParameter;
        this.stateDimension = stateDimension;
        this.stationaryMean = new double[stateDimension];
        this.compressedDGradient = new double[blockParameter.getCompressedDDimension()];
        this.rotationGradientFlat = new double[stateDimension * stateDimension];
        this.rotationGradient2D = new double[stateDimension][stateDimension];
        this.nativeGradient = new double[blockParameter.getBlockDiagonalNParameters()];
        this.diffusionGradient = new double[stateDimension * stateDimension];
        this.meanGradient = new double[stateDimension];
        this.smoothedInitialMean = new double[stateDimension];
        this.stateDiff = new double[stateDimension];
        this.initialCovInv = new double[stateDimension][stateDimension];
        this.smoothedInitialCovariance = new double[stateDimension][stateDimension];
        this.initialGradient = new double[stateDimension];
        this.preparedBasis = orthogonalParameterization.createPreparedBranchHandle();
        this.branchWorkspace = orthogonalParameterization.createBranchWorkspace();
        this.known = false;
    }

    public static boolean isAvailable(final OUProcessModel processModel,
                                      final Parameter selectionMatrixParameter) {
        return processModel != null
                && selectionMatrixParameter instanceof OrthogonalBlockDiagonalPolarStableMatrixParameter
                && processModel.getCovarianceGradientMethod()
                == OUProcessModel.CovarianceGradientMethod.STATIONARY_LYAPUNOV
                && processModel.getSelectionMatrixParameterization()
                instanceof OrthogonalBlockCanonicalParameterization;
    }

    public void makeDirty() {
        known = false;
    }

    public boolean supportsNativeSelectionParameter(final Parameter parameter) {
        if (parameter == blockParameter) {
            return false;
        }
        if (parameter == blockParameter.getParameter()
                || parameter == blockParameter.getRotationMatrixParameter()
                || parameter == blockParameter.getScalarBlockParameter()) {
            return true;
        }
        if (blockParameter.getRotationMatrixParameter() instanceof OrthogonalMatrixProvider
                && parameter == ((OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter())
                .getOrthogonalParameter()) {
            return true;
        }
        for (int i = 0; i < blockParameter.getTwoByTwoParameterFamilyCount(); ++i) {
            if (parameter == blockParameter.getTwoByTwoBlockParameter(i)) {
                return true;
            }
        }
        return false;
    }

    public boolean supportsDiffusionParameter(final Parameter parameter) {
        return diffusionParameterization.supportsParameter(parameter);
    }

    public boolean supportsMeanParameter(final Parameter parameter) {
        return parameter == stationaryMeanParameter;
    }

    public double[] getNativeSelectionGradient(final Parameter parameter,
                                               final CanonicalForwardTrajectory trajectory,
                                               final CanonicalBranchGradientCache branchCache,
                                               final TimeGrid timeGrid) {
        ensure(trajectory, branchCache, timeGrid);
        GaussianMatrixOps.copyFlatToMatrix(rotationGradientFlat, rotationGradient2D, stateDimension);
        blockParameter.chainGradient(compressedDGradient, nativeGradient);
        return assembleBlockGradientResult(parameter, nativeGradient, rotationGradient2D);
    }

    public double[] getDiffusionGradient(final Parameter parameter,
                                         final CanonicalForwardTrajectory trajectory,
                                         final CanonicalBranchGradientCache branchCache,
                                         final TimeGrid timeGrid) {
        ensure(trajectory, branchCache, timeGrid);
        if (parameter == diffusionParameterization.getMatrixParameter()) {
            return diffusionGradient.clone();
        }
        return diffusionParameterization.pullBackGradient(parameter, diffusionGradient);
    }

    public double[] getMeanGradient(final CanonicalForwardTrajectory trajectory,
                                    final CanonicalBranchGradientCache branchCache,
                                    final TimeGrid timeGrid) {
        ensure(trajectory, branchCache, timeGrid);
        if (stationaryMeanParameter.getDimension() == stateDimension) {
            return meanGradient.clone();
        }
        double sum = 0.0;
        for (int i = 0; i < stateDimension; ++i) {
            sum += meanGradient[i];
        }
        return new double[]{sum};
    }

    private void ensure(final CanonicalForwardTrajectory trajectory,
                        final CanonicalBranchGradientCache branchCache,
                        final TimeGrid timeGrid) {
        if (known) {
            return;
        }
        zero(compressedDGradient);
        zero(rotationGradientFlat);
        zero(diffusionGradient);
        zero(meanGradient);
        processModel.getInitialMean(stationaryMean);
        branchCache.ensure(trajectory);

        final MatrixParameterInterface diffusionMatrix = processModel.getDiffusionMatrix();
        for (int t = 0; t < trajectory.timeCount - 1; ++t) {
            final double dt = timeGrid.getDelta(t, t + 1);
            final CanonicalLocalTransitionAdjoints adjoints = branchCache.getAdjoints(t);
            orthogonalParameterization.prepareBranch(dt, stationaryMean, preparedBasis);
            orthogonalParameterization.accumulateNativeSelectionAndDiffusionGradientFromAdjointsPreparedFlat(
                    preparedBasis,
                    diffusionMatrix,
                    adjoints,
                    branchWorkspace,
                    compressedDGradient,
                    rotationGradientFlat,
                    false,
                    diffusionGradient);
            orthogonalParameterization.accumulateMeanGradientPrepared(
                    preparedBasis,
                    adjoints.dLogL_df,
                    meanGradient,
                    branchWorkspace);
        }
        accumulateInitialMeanGradient(trajectory);
        known = true;
    }

    private void accumulateInitialMeanGradient(final CanonicalForwardTrajectory trajectory) {
        GaussianFormConverter.fillMomentsFromState(
                trajectory.smoothedStates[0],
                smoothedInitialMean,
                smoothedInitialCovariance);
        GaussianMatrixOps.copyMatrix(initialCovarianceParameter.getParameterAsMatrix(), initialCovInv);
        final GaussianMatrixOps.CholeskyFactor initialChol =
                GaussianMatrixOps.cholesky(initialCovInv);
        GaussianMatrixOps.invertPositiveDefiniteFromCholesky(initialCovInv, initialChol);
        if (stationaryMeanParameter.getDimension() == 1) {
            final double meanValue = stationaryMeanParameter.getParameterValue(0);
            for (int i = 0; i < stateDimension; ++i) {
                stateDiff[i] = smoothedInitialMean[i] - meanValue;
            }
        } else {
            for (int i = 0; i < stateDimension; ++i) {
                stateDiff[i] = smoothedInitialMean[i] - stationaryMeanParameter.getParameterValue(i);
            }
        }
        GaussianMatrixOps.multiplyMatrixVector(initialCovInv, stateDiff, initialGradient);
        for (int i = 0; i < stateDimension; ++i) {
            meanGradient[i] += initialGradient[i];
        }
    }

    private double[] assembleBlockGradientResult(final Parameter requestedParameter,
                                                 final double[] nativeGradient,
                                                 final double[][] gradientR) {
        final int d = stateDimension;
        if (requestedParameter == blockParameter.getParameter()) {
            if (blockParameter.getRotationMatrixParameter() instanceof OrthogonalMatrixProvider) {
                final double[] angleGradient =
                        ((OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter())
                                .pullBackGradient(gradientR);
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
                && requestedParameter == ((OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter())
                .getOrthogonalParameter()) {
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

    private static void zero(final double[] vector) {
        for (int i = 0; i < vector.length; ++i) {
            vector[i] = 0.0;
        }
    }

    private static double[] flattenColumnMajor(final double[][] matrix) {
        final double[] out = new double[matrix.length * matrix.length];
        flattenColumnMajor(matrix, out, 0);
        return out;
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
