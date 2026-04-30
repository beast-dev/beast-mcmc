package dr.inference.timeseries.engine.gaussian;

import dr.evomodel.treedatalikelihood.continuous.gaussian.message.GaussianMatrixOps;

import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalBranchMessageContributionUtils;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalTransitionAdjointUtils;

import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.OrthogonalBlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.OrthogonalMatrixProvider;
import dr.inference.model.Parameter;
import dr.inference.timeseries.core.TimeGrid;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockCanonicalParameterization;
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

    private final CanonicalBranchMessageContribution localContribution;
    private final CanonicalBranchMessageContributionUtils.Workspace contributionWorkspace;
    private final CanonicalLocalTransitionAdjoints localAdjoints;
    private final CanonicalTransitionAdjointUtils.Workspace canonicalAdjointWorkspace;
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

        this.localContribution = new CanonicalBranchMessageContribution(stateDimension);
        this.contributionWorkspace = new CanonicalBranchMessageContributionUtils.Workspace(stateDimension);
        this.localAdjoints = new CanonicalLocalTransitionAdjoints(stateDimension);
        this.canonicalAdjointWorkspace = new CanonicalTransitionAdjointUtils.Workspace(stateDimension);
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
            CanonicalBranchMessageContributionUtils.fillFromPairState(
                    trajectory.branchPairStates[t],
                    contributionWorkspace,
                    localContribution);
            CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                    trajectory.transitions[t],
                    localContribution,
                    canonicalAdjointWorkspace,
                    localAdjoints);
            GaussianMatrixOps.copyFlatToMatrix(localAdjoints.dLogL_dF, dLogL_dF, d);
            repr.accumulateSelectionGradient(
                    t, t + 1, timeGrid, dLogL_dF, localAdjoints.dLogL_df, gradientAccumulator);

            GaussianMatrixOps.copyFlatToMatrix(localAdjoints.dLogL_dOmega, dLogL_dV, d);
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

    private double[] computeOrthogonalNativeGradient(final Parameter requestedParameter,
                                                     final CanonicalForwardTrajectory trajectory,
                                                     final TimeGrid timeGrid) {
        final OrthogonalBlockCanonicalParameterization orthogonalParameterization =
                (OrthogonalBlockCanonicalParameterization) processModel.getSelectionMatrixParameterization();
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
                    contributionWorkspace,
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
