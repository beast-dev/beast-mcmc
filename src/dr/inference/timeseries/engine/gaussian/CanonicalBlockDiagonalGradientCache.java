package dr.inference.timeseries.engine.gaussian;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.DiffusionMatrixParameterization;
import dr.evomodel.continuous.ou.canonical.CanonicalBranchWorkspace;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedBranchGradientAccumulator;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedBranchHandle;
import dr.evomodel.continuous.ou.blockdiagonal.BlockDiagonalNativeCanonicalParameterization;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.GaussianFormConverter;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;
import dr.inference.timeseries.representation.KernelBackedGaussianTransitionRepresentation;

/**
 * Timeseries-only joint accumulator for block-diagonal canonical OU gradients.
 */
public final class CanonicalBlockDiagonalGradientCache {

    private final OUProcessModel processModel;
    private final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter;
    private final BlockDiagonalNativeCanonicalParameterization blockParameterization;
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
    private final double[] dBasisMeanGradient;
    private final double[] smoothedInitialMean;
    private final double[] stateDiff;
    private final double[] flatSmoothedInitialCovariance;
    private final double[] flatInitialCovariance;
    private final double[] flatInitialCovarianceInverse;
    private final double[] flatInitialCholesky;
    private final double[] flatInitialLowerInverse;
    private final double[] initialGradient;
    private final GaussianFormConverter.Workspace converterWorkspace;
    private final CanonicalPreparedBranchGradientAccumulator preparedGradientAccumulator;
    private final CanonicalPreparedBranchHandle preparedBasis;
    private final CanonicalBranchWorkspace branchWorkspace;
    private boolean selectionAndMeanKnown;
    private boolean diffusionKnown;

    public CanonicalBlockDiagonalGradientCache(final OUProcessModel processModel,
                                               final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                               final DiffusionMatrixParameterization diffusionParameterization,
                                               final Parameter stationaryMeanParameter,
                                               final MatrixParameter initialCovarianceParameter,
                                               final int stateDimension) {
        this.processModel = processModel;
        this.blockParameter = blockParameter;
        this.blockParameterization =
                (BlockDiagonalNativeCanonicalParameterization) processModel.getSelectionMatrixParameterization();
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
        this.dBasisMeanGradient = new double[stateDimension];
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
        this.preparedGradientAccumulator =
                new CanonicalPreparedBranchGradientAccumulator(blockParameterization, stateDimension);
        this.preparedBasis = blockParameterization.createPreparedBranchHandle();
        this.branchWorkspace = blockParameterization.createBranchWorkspace();
        this.preparedGradientAccumulator.checkBuffers(
                compressedDGradient,
                nativeGradient,
                rotationGradientFlat,
                dBasisDiffusionGradient,
                dBasisMeanGradient);
        this.selectionAndMeanKnown = false;
        this.diffusionKnown = false;
    }

    public static boolean isAvailable(final OUProcessModel processModel,
                                      final Parameter selectionMatrixParameter) {
        return BlockDiagonalFormulaSupport.isNativeSelectionAvailable(
                processModel,
                selectionMatrixParameter);
    }

    public void makeDirty() {
        selectionAndMeanKnown = false;
        diffusionKnown = false;
    }

    public boolean supportsNativeSelectionParameter(final Parameter parameter) {
        return BlockDiagonalFormulaSupport.supportsNativeSelectionParameter(
                blockParameter,
                parameter);
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
        BlockDiagonalFormulaSupport.zero(compressedDGradient);
        BlockDiagonalFormulaSupport.zero(rotationGradientFlat);
        BlockDiagonalFormulaSupport.zero(meanGradient);
        BlockDiagonalFormulaSupport.zero(dBasisMeanGradient);
        processModel.getInitialMean(stationaryMean);
        branchCache.ensure(trajectory);

        final MatrixParameterInterface diffusionMatrix = processModel.getDiffusionMatrix();
        for (int t = 0; t < trajectory.timeCount - 1; ++t) {
            final double dt = timeGrid.getDelta(t, t + 1);
            final CanonicalLocalTransitionAdjoints adjoints = branchCache.getAdjoints(t);
            final CanonicalPreparedBranchHandle prepared =
                    threadPreparedBranch(transitionRepresentation, dt, stationaryMean);
            preparedGradientAccumulator.accumulateSelectionAndMean(
                    prepared,
                    diffusionMatrix,
                    adjoints,
                    branchWorkspace,
                    compressedDGradient,
                    rotationGradientFlat,
                    meanGradient,
                    dBasisMeanGradient);
        }
        preparedGradientAccumulator.finishDelayedMean(
                dBasisMeanGradient,
                meanGradient,
                branchWorkspace);
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
                && blockParameterization.supportsDelayedDiffusionGradientRotation()) {
            ensureSelectionMeanAndDiffusion(trajectory, branchCache, transitionRepresentation, timeGrid);
            return;
        }
        BlockDiagonalFormulaSupport.zero(diffusionGradient);
        BlockDiagonalFormulaSupport.zero(dBasisDiffusionGradient);
        processModel.getInitialMean(stationaryMean);
        branchCache.ensure(trajectory);

        for (int t = 0; t < trajectory.timeCount - 1; ++t) {
            final double dt = timeGrid.getDelta(t, t + 1);
            final CanonicalLocalTransitionAdjoints adjoints = branchCache.getAdjoints(t);
            final CanonicalPreparedBranchHandle prepared =
                    threadPreparedBranch(transitionRepresentation, dt, stationaryMean);
            preparedGradientAccumulator.accumulateDiffusion(
                    prepared,
                    adjoints,
                    branchWorkspace,
                    false,
                    false,
                    diffusionGradient,
                    dBasisDiffusionGradient);
        }
        preparedGradientAccumulator.finishDelayedDiffusion(
                dBasisDiffusionGradient,
                diffusionGradient,
                branchWorkspace);
        diffusionKnown = true;
    }

    private void ensureSelectionMeanAndDiffusion(final CanonicalForwardTrajectory trajectory,
                                                 final CanonicalBranchGradientCache branchCache,
                                                 final GaussianTransitionRepresentation transitionRepresentation,
                                                 final TimeGrid timeGrid) {
        BlockDiagonalFormulaSupport.zero(compressedDGradient);
        BlockDiagonalFormulaSupport.zero(rotationGradientFlat);
        BlockDiagonalFormulaSupport.zero(meanGradient);
        BlockDiagonalFormulaSupport.zero(diffusionGradient);
        BlockDiagonalFormulaSupport.zero(dBasisDiffusionGradient);
        BlockDiagonalFormulaSupport.zero(dBasisMeanGradient);
        processModel.getInitialMean(stationaryMean);
        branchCache.ensure(trajectory);

        final MatrixParameterInterface diffusionMatrix = processModel.getDiffusionMatrix();
        for (int t = 0; t < trajectory.timeCount - 1; ++t) {
            final double dt = timeGrid.getDelta(t, t + 1);
            final CanonicalLocalTransitionAdjoints adjoints = branchCache.getAdjoints(t);
            final CanonicalPreparedBranchHandle prepared =
                    threadPreparedBranch(transitionRepresentation, dt, stationaryMean);
            preparedGradientAccumulator.accumulateSelectionDiffusionAndMean(
                    prepared,
                    diffusionMatrix,
                    adjoints,
                    branchWorkspace,
                    compressedDGradient,
                    rotationGradientFlat,
                    diffusionGradient,
                    dBasisDiffusionGradient,
                    meanGradient,
                    dBasisMeanGradient);
        }
        preparedGradientAccumulator.finishDelayedDiffusion(
                dBasisDiffusionGradient,
                diffusionGradient,
                branchWorkspace);
        preparedGradientAccumulator.finishDelayedMean(
                dBasisMeanGradient,
                meanGradient,
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
            final CanonicalPreparedBranchHandle cached =
                    ((KernelBackedGaussianTransitionRepresentation) transitionRepresentation)
                            .getReusablePreparedCanonicalBranch(dt);
            if (cached != null) {
                return cached;
            }
            final CanonicalPreparedBranchHandle prepared =
                    ((KernelBackedGaussianTransitionRepresentation) transitionRepresentation)
                            .getThreadPreparedCanonicalBranch(dt, stationaryMean);
            if (prepared != null) {
                return prepared;
            }
        }
        blockParameterization.prepareBranch(dt, stationaryMean, preparedBasis);
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
        return BlockDiagonalFormulaSupport.assembleNativeSelectionGradient(
                stateDimension,
                requestedParameter,
                blockParameter,
                nativeGradient,
                flatGradientR);
    }
}
