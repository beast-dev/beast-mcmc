package dr.inference.timeseries.engine.gaussian;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianMatrixOps;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContributionUtils;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalTransitionAdjointUtils;

import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.OrthogonalBlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.OrthogonalMatrixProvider;
import dr.inference.model.Parameter;
import dr.inference.timeseries.core.TimeGrid;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockCanonicalParameterization;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;

import java.util.Arrays;

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

    private final CanonicalBranchMessageContribution localContribution;
    private final CanonicalBranchMessageContributionUtils.Workspace contributionWorkspace;
    private final CanonicalLocalTransitionAdjoints localAdjoints;
    private final CanonicalTransitionAdjointUtils.Workspace canonicalAdjointWorkspace;
    private final double[][] dLogL_dF;
    private final double[][] dLogL_dV;
    private final double[] nativeCompressedGradientScratch;
    // Row-major flat rotation gradient (replaces former double[][] field).
    private final double[] nativeRotationGradientFlat;
    // 2-D scratch for pullBackGradient() calls that still require double[][].
    private final double[][] nativeRotationGradient2D;
    private final double[] stationaryMeanScratch;
    private final CanonicalOrthogonalBlockGradientCache orthogonalBlockGradientCache;

    public CanonicalSelectionMatrixGradientFormula(final Parameter selectionMatrixParameter,
                                                   final int stateDimension) {
        this(null, selectionMatrixParameter, stateDimension);
    }

    public CanonicalSelectionMatrixGradientFormula(final OUProcessModel processModel,
                                                   final Parameter selectionMatrixParameter,
                                                   final int stateDimension) {
        this(processModel, selectionMatrixParameter, stateDimension, null);
    }

    public CanonicalSelectionMatrixGradientFormula(final OUProcessModel processModel,
                                                   final Parameter selectionMatrixParameter,
                                                   final int stateDimension,
                                                   final CanonicalOrthogonalBlockGradientCache orthogonalBlockGradientCache) {
        if (selectionMatrixParameter == null) {
            throw new IllegalArgumentException("selectionMatrixParameter must not be null");
        }
        if (stateDimension < 1) {
            throw new IllegalArgumentException("stateDimension must be at least 1");
        }
        this.selectionMatrixParameter = selectionMatrixParameter;
        this.stateDimension = stateDimension;
        this.processModel = processModel;
        this.orthogonalBlockGradientCache = orthogonalBlockGradientCache;

        this.localContribution = new CanonicalBranchMessageContribution(stateDimension);
        this.contributionWorkspace = new CanonicalBranchMessageContributionUtils.Workspace(stateDimension);
        this.localAdjoints = new CanonicalLocalTransitionAdjoints(stateDimension);
        this.canonicalAdjointWorkspace = new CanonicalTransitionAdjointUtils.Workspace(stateDimension);
        this.dLogL_dF = new double[stateDimension][stateDimension];
        this.dLogL_dV = new double[stateDimension][stateDimension];
        final boolean hasBlock = selectionMatrixParameter instanceof AbstractBlockDiagonalTwoByTwoMatrixParameter;
        this.nativeCompressedGradientScratch = hasBlock
                ? new double[((AbstractBlockDiagonalTwoByTwoMatrixParameter) selectionMatrixParameter).getCompressedDDimension()]
                : null;
        this.nativeRotationGradientFlat = hasBlock
                ? new double[stateDimension * stateDimension]
                : null;
        this.nativeRotationGradient2D = hasBlock
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
        return computeGradient(parameter, trajectory, null, repr, timeGrid);
    }

    @Override
    public double[] computeGradient(final Parameter parameter,
                                    final CanonicalForwardTrajectory trajectory,
                                    final CanonicalBranchGradientCache branchGradientCache,
                                    final GaussianTransitionRepresentation repr,
                                    final TimeGrid timeGrid) {
        if (shouldUseOrthogonalNativePath(parameter)) {
            if (orthogonalBlockGradientCache != null
                    && orthogonalBlockGradientCache.supportsNativeSelectionParameter(parameter)
                    && branchGradientCache != null) {
                return orthogonalBlockGradientCache.getNativeSelectionGradient(
                        parameter, trajectory, branchGradientCache, timeGrid);
            }
            return computeOrthogonalNativeGradient(parameter, trajectory, branchGradientCache, timeGrid);
        }

        final int d = stateDimension;
        final int T = trajectory.timeCount;
        final double[] gradientAccumulator = new double[d * d];
        if (branchGradientCache != null) {
            branchGradientCache.ensure(trajectory);
        }

        for (int t = 0; t < T - 1; ++t) {
            final CanonicalLocalTransitionAdjoints adjoints =
                    localAdjoints(t, trajectory, branchGradientCache);
            GaussianMatrixOps.copyFlatToMatrix(adjoints.dLogL_dF, dLogL_dF, d);
            repr.accumulateSelectionGradient(
                    t, t + 1, timeGrid, dLogL_dF, adjoints.dLogL_df, gradientAccumulator);

            GaussianMatrixOps.copyFlatToMatrix(adjoints.dLogL_dOmega, dLogL_dV, d);
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
                instanceof OrthogonalBlockCanonicalParameterization;
    }

    @Override
    public void makeDirty() {
        if (orthogonalBlockGradientCache != null) {
            orthogonalBlockGradientCache.makeDirty();
        }
    }

    private double[] computeOrthogonalNativeGradient(final Parameter requestedParameter,
                                                     final CanonicalForwardTrajectory trajectory,
                                                     final CanonicalBranchGradientCache branchGradientCache,
                                                     final TimeGrid timeGrid) {
        final OrthogonalBlockCanonicalParameterization orthogonalParameterization =
                (OrthogonalBlockCanonicalParameterization) processModel.getSelectionMatrixParameterization();
        final OrthogonalBlockDiagonalPolarStableMatrixParameter blockParameter =
                (OrthogonalBlockDiagonalPolarStableMatrixParameter) selectionMatrixParameter;
        final int T = trajectory.timeCount;
        final double[] compressedDGradient = nativeCompressedGradientScratch;
        final double[] rotationGradientFlat = nativeRotationGradientFlat;
        final double[] stationaryMean = stationaryMeanScratch;
        Arrays.fill(compressedDGradient, 0.0);
        Arrays.fill(rotationGradientFlat, 0.0);
        processModel.getInitialMean(stationaryMean);
        if (branchGradientCache != null) {
            branchGradientCache.ensure(trajectory);
        }

        for (int t = 0; t < T - 1; ++t) {
            final double dt = timeGrid.getDelta(t, t + 1);
            final CanonicalBranchMessageContribution contribution =
                    localContribution(t, trajectory, branchGradientCache);
            orthogonalParameterization.accumulateNativeGradientFromCanonicalContributionFlat(
                    processModel.getDiffusionMatrix(),
                    stationaryMean,
                    dt,
                    contribution,
                    localAdjoints(t, trajectory, branchGradientCache),
                    compressedDGradient,
                    rotationGradientFlat);
        }

        // Convert flat rotation gradient to 2D for assembleBlockGradientResult.
        GaussianMatrixOps.copyFlatToMatrix(rotationGradientFlat, nativeRotationGradient2D, stateDimension);

        final double[] nativeGradient = new double[blockParameter.getBlockDiagonalNParameters()];
        blockParameter.chainGradient(compressedDGradient, nativeGradient);
        return assembleBlockGradientResult(requestedParameter, blockParameter, nativeGradient, nativeRotationGradient2D);
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

    private CanonicalBranchMessageContribution localContribution(final int branchIndex,
                                                                 final CanonicalForwardTrajectory trajectory,
                                                                 final CanonicalBranchGradientCache branchGradientCache) {
        if (branchGradientCache != null) {
            return branchGradientCache.getContribution(branchIndex);
        }
        CanonicalBranchMessageContributionUtils.fillFromPairState(
                trajectory.branchPairStates[branchIndex],
                contributionWorkspace,
                localContribution);
        return localContribution;
    }

    private CanonicalLocalTransitionAdjoints localAdjoints(final int branchIndex,
                                                           final CanonicalForwardTrajectory trajectory,
                                                           final CanonicalBranchGradientCache branchGradientCache) {
        if (branchGradientCache != null) {
            return branchGradientCache.getAdjoints(branchIndex);
        }
        CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                trajectory.transitions[branchIndex],
                localContribution(branchIndex, trajectory, null),
                canonicalAdjointWorkspace,
                localAdjoints);
        return localAdjoints;
    }
}
