package dr.inference.timeseries.engine.kalman.formula;

import dr.evomodel.continuous.ou.DiffusionMatrixParameterization;
import dr.inference.model.Parameter;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.engine.kalman.BranchSmootherStats;
import dr.inference.timeseries.engine.kalman.ForwardTrajectory;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;

/**
 * Analytical gradient of the Kalman log-likelihood with respect to the diffusion
 * matrix Q (or native parameters that determine Q).
 */
public final class ExpectationDiffusionMatrixGradientFormula implements ExpectationGradientFormula {

    private final DiffusionMatrixParameterization diffusionParameterization;
    private final int stateDimension;
    private final ExpectationGaussianBranchGradientAdjoints branchAdjoints;

    public ExpectationDiffusionMatrixGradientFormula(final DiffusionMatrixParameterization diffusionParameterization,
                                          final int stateDimension) {
        if (diffusionParameterization == null) {
            throw new IllegalArgumentException("diffusionParameterization must not be null");
        }
        if (stateDimension < 1) {
            throw new IllegalArgumentException("stateDimension must be at least 1");
        }
        this.diffusionParameterization = diffusionParameterization;
        this.stateDimension = stateDimension;
        this.branchAdjoints = new ExpectationGaussianBranchGradientAdjoints(stateDimension);
    }

    @Override
    public boolean supportsParameter(final Parameter parameter) {
        return diffusionParameterization.supportsParameter(parameter);
    }

    @Override
    public double[] computeGradient(final Parameter parameter,
                                    final BranchSmootherStats[] smootherStats,
                                    final ForwardTrajectory trajectory,
                                    final GaussianTransitionRepresentation repr,
                                    final TimeGrid timeGrid) {
        final int T = trajectory.timeCount;
        final int d = stateDimension;
        final double[] gradientAccumulator = new double[d * d];

        for (int t = 0; t < T - 1; ++t) {
            final int transitionMatrixOffset = trajectory.branchMatrixOffset(t);
            final int transitionOffsetOffset = trajectory.branchVectorOffset(t);
            branchAdjoints.compute(
                    smootherStats[t],
                    smootherStats[t + 1],
                    trajectory.transitionMatrices,
                    transitionMatrixOffset,
                    trajectory.transitionOffsets,
                    transitionOffsetOffset,
                    trajectory.stepCovariances,
                    transitionMatrixOffset);

            repr.accumulateDiffusionGradientFlat(t, t + 1, timeGrid,
                    branchAdjoints.dLogL_dV(), gradientAccumulator);
        }

        return diffusionParameterization.pullBackGradient(parameter, gradientAccumulator);
    }
}
