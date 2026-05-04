package dr.inference.timeseries.engine.gaussian;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.canonical.CanonicalBranchWorkspace;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedBranchHandle;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockCanonicalParameterization;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.GaussianFormConverter;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.OrthogonalBlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.OrthogonalMatrixProvider;
import dr.inference.model.Parameter;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.gaussian.DiffusionMatrixParameterization;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;
import dr.inference.timeseries.representation.KernelBackedGaussianTransitionRepresentation;

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
    private final double[] nativeGradient;
    private final double[] diffusionGradient;
    private final double[] dBasisDiffusionGradient;
    private final double[] meanGradient;
    private final double[] smoothedInitialMean;
    private final double[] stateDiff;
    private final double[] flatSmoothedInitialCovariance;
    private final double[] flatInitialCovariance;
    private final double[] flatInitialCovarianceInverse;
    private final double[] flatInitialCholesky;
    private final double[] flatInitialLowerInverse;
    private final double[] initialGradient;
    private final GaussianFormConverter.Workspace converterWorkspace;
    private final CanonicalPreparedBranchHandle preparedBasis;
    private final CanonicalBranchWorkspace branchWorkspace;
    private boolean selectionAndMeanKnown;
    private boolean diffusionKnown;

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
        this.nativeGradient = new double[blockParameter.getBlockDiagonalNParameters()];
        this.diffusionGradient = new double[stateDimension * stateDimension];
        this.dBasisDiffusionGradient = new double[stateDimension * stateDimension];
        this.meanGradient = new double[stateDimension];
        this.smoothedInitialMean = new double[stateDimension];
        this.stateDiff = new double[stateDimension];
        this.flatSmoothedInitialCovariance = new double[stateDimension * stateDimension];
        this.flatInitialCovariance = new double[stateDimension * stateDimension];
        this.flatInitialCovarianceInverse = new double[stateDimension * stateDimension];
        this.flatInitialCholesky = new double[stateDimension * stateDimension];
        this.flatInitialLowerInverse = new double[stateDimension * stateDimension];
        this.initialGradient = new double[stateDimension];
        this.converterWorkspace = new GaussianFormConverter.Workspace();
        this.converterWorkspace.ensureDim(stateDimension);
        this.preparedBasis = orthogonalParameterization.createPreparedBranchHandle();
        this.branchWorkspace = orthogonalParameterization.createBranchWorkspace();
        this.selectionAndMeanKnown = false;
        this.diffusionKnown = false;
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
        selectionAndMeanKnown = false;
        diffusionKnown = false;
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
                                               final GaussianTransitionRepresentation transitionRepresentation,
                                               final TimeGrid timeGrid) {
        ensureSelectionAndMean(trajectory, branchCache, transitionRepresentation, timeGrid);
        blockParameter.chainGradient(compressedDGradient, nativeGradient);
        return assembleBlockGradientResult(parameter, nativeGradient, rotationGradientFlat);
    }

    public double[] getDiffusionGradient(final Parameter parameter,
                                         final CanonicalForwardTrajectory trajectory,
                                         final CanonicalBranchGradientCache branchCache,
                                         final GaussianTransitionRepresentation transitionRepresentation,
                                         final TimeGrid timeGrid) {
        ensureDiffusion(trajectory, branchCache, transitionRepresentation, timeGrid);
        if (parameter == diffusionParameterization.getMatrixParameter()) {
            return diffusionGradient.clone();
        }
        return diffusionParameterization.pullBackGradient(parameter, diffusionGradient);
    }

    public double[] getMeanGradient(final CanonicalForwardTrajectory trajectory,
                                    final CanonicalBranchGradientCache branchCache,
                                    final GaussianTransitionRepresentation transitionRepresentation,
                                    final TimeGrid timeGrid) {
        ensureSelectionAndMean(trajectory, branchCache, transitionRepresentation, timeGrid);
        if (stationaryMeanParameter.getDimension() == stateDimension) {
            return meanGradient.clone();
        }
        double sum = 0.0;
        for (int i = 0; i < stateDimension; ++i) {
            sum += meanGradient[i];
        }
        return new double[]{sum};
    }

    private void ensureSelectionAndMean(final CanonicalForwardTrajectory trajectory,
                                        final CanonicalBranchGradientCache branchCache,
                                        final GaussianTransitionRepresentation transitionRepresentation,
                                        final TimeGrid timeGrid) {
        if (selectionAndMeanKnown) {
            return;
        }
        zero(compressedDGradient);
        zero(rotationGradientFlat);
        zero(meanGradient);
        processModel.getInitialMean(stationaryMean);
        branchCache.ensure(trajectory);

        final MatrixParameterInterface diffusionMatrix = processModel.getDiffusionMatrix();
        for (int t = 0; t < trajectory.timeCount - 1; ++t) {
            final double dt = timeGrid.getDelta(t, t + 1);
            final CanonicalLocalTransitionAdjoints adjoints = branchCache.getAdjoints(t);
            final CanonicalPreparedBranchHandle prepared =
                    threadPreparedBranch(transitionRepresentation, dt, stationaryMean);
            orthogonalParameterization.accumulateNativeGradientFromAdjointsPreparedFlat(
                    prepared,
                    diffusionMatrix,
                    adjoints,
                    branchWorkspace,
                    compressedDGradient,
                    rotationGradientFlat);
            orthogonalParameterization.accumulateMeanGradientPrepared(
                    prepared,
                    adjoints.dLogL_df,
                    meanGradient,
                    branchWorkspace);
        }
        accumulateInitialMeanGradient(trajectory);
        selectionAndMeanKnown = true;
    }

    private void ensureDiffusion(final CanonicalForwardTrajectory trajectory,
                                 final CanonicalBranchGradientCache branchCache,
                                 final GaussianTransitionRepresentation transitionRepresentation,
                                 final TimeGrid timeGrid) {
        if (diffusionKnown) {
            return;
        }
        if (!selectionAndMeanKnown
                && orthogonalParameterization.supportsDelayedDiffusionGradientRotation()) {
            ensureSelectionMeanAndDiffusion(trajectory, branchCache, transitionRepresentation, timeGrid);
            return;
        }
        zero(diffusionGradient);
        processModel.getInitialMean(stationaryMean);
        branchCache.ensure(trajectory);

        for (int t = 0; t < trajectory.timeCount - 1; ++t) {
            final double dt = timeGrid.getDelta(t, t + 1);
            final CanonicalLocalTransitionAdjoints adjoints = branchCache.getAdjoints(t);
            final CanonicalPreparedBranchHandle prepared =
                    threadPreparedBranch(transitionRepresentation, dt, stationaryMean);
            orthogonalParameterization.accumulateDiffusionGradientPreparedFlat(
                    prepared,
                    adjoints.dLogL_dOmega,
                    false,
                    diffusionGradient,
                    branchWorkspace);
        }
        diffusionKnown = true;
    }

    private void ensureSelectionMeanAndDiffusion(final CanonicalForwardTrajectory trajectory,
                                                 final CanonicalBranchGradientCache branchCache,
                                                 final GaussianTransitionRepresentation transitionRepresentation,
                                                 final TimeGrid timeGrid) {
        zero(compressedDGradient);
        zero(rotationGradientFlat);
        zero(meanGradient);
        zero(diffusionGradient);
        zero(dBasisDiffusionGradient);
        processModel.getInitialMean(stationaryMean);
        branchCache.ensure(trajectory);

        final MatrixParameterInterface diffusionMatrix = processModel.getDiffusionMatrix();
        for (int t = 0; t < trajectory.timeCount - 1; ++t) {
            final double dt = timeGrid.getDelta(t, t + 1);
            final CanonicalLocalTransitionAdjoints adjoints = branchCache.getAdjoints(t);
            final CanonicalPreparedBranchHandle prepared =
                    threadPreparedBranch(transitionRepresentation, dt, stationaryMean);
            orthogonalParameterization.accumulateNativeSelectionAndDiffusionGradientFromAdjointsPreparedFlat(
                    prepared,
                    diffusionMatrix,
                    adjoints,
                    branchWorkspace,
                    compressedDGradient,
                    rotationGradientFlat,
                    true,
                    dBasisDiffusionGradient);
            orthogonalParameterization.accumulateMeanGradientPrepared(
                    prepared,
                    adjoints.dLogL_df,
                    meanGradient,
                    branchWorkspace);
        }
        orthogonalParameterization.finishDiffusionGradientFromDBasisFlat(
                dBasisDiffusionGradient,
                diffusionGradient,
                branchWorkspace);
        accumulateInitialMeanGradient(trajectory);
        selectionAndMeanKnown = true;
        diffusionKnown = true;
    }

    private CanonicalPreparedBranchHandle threadPreparedBranch(
            final GaussianTransitionRepresentation transitionRepresentation,
            final double dt,
            final double[] stationaryMean) {
        if (transitionRepresentation instanceof KernelBackedGaussianTransitionRepresentation) {
            final CanonicalPreparedBranchHandle prepared =
                    ((KernelBackedGaussianTransitionRepresentation) transitionRepresentation)
                            .getThreadPreparedCanonicalBranch(dt, stationaryMean);
            if (prepared != null) {
                return prepared;
            }
        }
        orthogonalParameterization.prepareBranch(dt, stationaryMean, preparedBasis);
        return preparedBasis;
    }

    private void accumulateInitialMeanGradient(final CanonicalForwardTrajectory trajectory) {
        GaussianFormConverter.fillMomentsFromState(
                trajectory.smoothedStates[0],
                smoothedInitialMean,
                flatSmoothedInitialCovariance,
                stateDimension,
                converterWorkspace);
        fillInitialCovarianceFlat();
        if (!MatrixOps.tryCholesky(flatInitialCovariance, flatInitialCholesky, stateDimension)) {
            throw new IllegalArgumentException("Initial covariance matrix is not positive definite");
        }
        MatrixOps.invertFromCholesky(
                flatInitialCholesky,
                flatInitialLowerInverse,
                flatInitialCovarianceInverse,
                stateDimension);
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
        MatrixOps.matVec(flatInitialCovarianceInverse, stateDiff, initialGradient, stateDimension);
        for (int i = 0; i < stateDimension; ++i) {
            meanGradient[i] += initialGradient[i];
        }
    }

    private void fillInitialCovarianceFlat() {
        for (int row = 0; row < stateDimension; ++row) {
            final int rowOffset = row * stateDimension;
            for (int col = 0; col < stateDimension; ++col) {
                flatInitialCovariance[rowOffset + col] =
                        initialCovarianceParameter.getParameterValue(row, col);
            }
        }
    }

    private double[] assembleBlockGradientResult(final Parameter requestedParameter,
                                                 final double[] nativeGradient,
                                                 final double[] flatGradientR) {
        final int d = stateDimension;
        if (requestedParameter == blockParameter.getParameter()) {
            if (blockParameter.getRotationMatrixParameter() instanceof OrthogonalMatrixProvider) {
                final OrthogonalMatrixProvider provider =
                        (OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter();
                final int angleCount = provider.getOrthogonalParameter().getDimension();
                final double[] out = new double[nativeGradient.length + angleCount];
                System.arraycopy(nativeGradient, 0, out, 0, nativeGradient.length);
                provider.fillPullBackGradientFlat(flatGradientR, d, out, nativeGradient.length);
                return out;
            }
            final double[] out = new double[nativeGradient.length + d * d];
            System.arraycopy(nativeGradient, 0, out, 0, nativeGradient.length);
            flattenColumnMajor(flatGradientR, out, nativeGradient.length, d);
            return out;
        }
        if (requestedParameter == blockParameter.getRotationMatrixParameter()) {
            return flattenColumnMajor(flatGradientR, d);
        }
        if (blockParameter.getRotationMatrixParameter() instanceof OrthogonalMatrixProvider
                && requestedParameter == ((OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter())
                .getOrthogonalParameter()) {
            final double[] out = new double[requestedParameter.getDimension()];
            ((OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter())
                    .fillPullBackGradientFlat(flatGradientR, d, out);
            return out;
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

    private static double[] flattenColumnMajor(final double[] rowMajor,
                                               final int dimension) {
        final double[] out = new double[dimension * dimension];
        flattenColumnMajor(rowMajor, out, 0, dimension);
        return out;
    }

    private static void flattenColumnMajor(final double[] rowMajor,
                                           final double[] out,
                                           final int offset,
                                           final int dimension) {
        int index = offset;
        for (int col = 0; col < dimension; ++col) {
            for (int row = 0; row < dimension; ++row) {
                out[index++] = rowMajor[row * dimension + col];
            }
        }
    }
}
