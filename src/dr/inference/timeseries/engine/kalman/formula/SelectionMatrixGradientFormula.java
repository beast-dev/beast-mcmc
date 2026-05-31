package dr.inference.timeseries.engine.kalman.formula;

import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.Parameter;
import dr.evomodel.treedatalikelihood.continuous.canonical.gradient.CanonicalSelectionGradientProjector;
import dr.inference.timeseries.core.TimeGrid;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.blockdiagonal.BlockDiagonalNativeCanonicalParameterization;
import dr.inference.timeseries.engine.kalman.BranchSmootherStats;
import dr.inference.timeseries.engine.kalman.ForwardTrajectory;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;

import java.util.Arrays;

/**
 * Analytical gradient of log p(Y | θ) with respect to the selection (drift) matrix A.
 *
 * <h3>Mathematical derivation</h3>
 * By Fisher's identity the gradient decomposes as a branch sum:
 * <pre>
 *   ∂logL/∂A  =  Σ_t  E_{x|Y} [ ∂ log p(x_{t+1} | x_t; A) / ∂A ]
 * </pre>
 * At each branch t the expectation over the smoothed trajectory yields the
 * <em>transition-space sensitivities</em>
 * <pre>
 *   ∂logL/∂F_t  =  Q_step_t^{-1} · B_t
 *   ∂logL/∂f_t  =  Q_step_t^{-1} · r̄_t
 * </pre>
 * where
 * <pre>
 *   B_t  =  V_{t+1,t}  −  F_t · P_{t|T}  +  r̄_t · m_{t|T}^T
 *   r̄_t  =  m_{t+1|T}  −  F_t · m_{t|T}  −  f_t
 * </pre>
 * and the lag-1 cross-covariance is recovered without extra storage as
 * <pre>
 *   V_{t+1,t}  =  P_{t+1|T} · G_t^T
 * </pre>
 * The formula accumulates two contributions to ∂logL/∂A per branch:
 * <ol>
 *   <li><b>F/f-path</b>: via ∂logL/∂F_t and ∂logL/∂f_t, delegated to
 *       {@link GaussianTransitionRepresentation#accumulateSelectionGradient}.</li>
 *   <li><b>V-path</b>: via ∂logL/∂V_t computed from the expected squared prediction
 *       error S_t = E[δ_t δ_t^T | Y] and the step-covariance inverse, delegated to
 *       {@link GaussianTransitionRepresentation#accumulateSelectionGradientFromCovariance}.
 *       For models where V does not depend on A (Euler), this contribution is zero.</li>
 * </ol>
 *
 * <h3>Complexity</h3>
 * One forward + backward smoother pass suffices for all parameter gradients.  This
 * formula contributes O(T · d³) additional work on top of that shared pass.
 */
public class SelectionMatrixGradientFormula implements GradientFormula {

    private final Parameter selectionMatrixParameter;
    private final int stateDimension;
    private final CanonicalSelectionGradientProjector.Workspace blockProjectorWorkspace;

    private final GaussianBranchGradientAdjoints branchAdjoints;
    private final double[]   dLogL_dFFlat;
    private final double[]   dLogL_dVFlat;
    private final double[]   rotationGradientFlat;

    /**
     * @param selectionMatrixParameter the drift/selection matrix {@link Parameter} whose
     *                                 gradient is computed; compared by identity in
     *                                 {@link #supportsParameter}
     * @param stateDimension           state-space dimension d
     */
    public SelectionMatrixGradientFormula(final Parameter selectionMatrixParameter,
                                          final int stateDimension) {
        if (selectionMatrixParameter == null) {
            throw new IllegalArgumentException("selectionMatrixParameter must not be null");
        }
        if (stateDimension < 1) {
            throw new IllegalArgumentException("stateDimension must be at least 1");
        }
        this.selectionMatrixParameter = selectionMatrixParameter;
        this.stateDimension           = stateDimension;
        this.blockProjectorWorkspace = selectionMatrixParameter instanceof AbstractBlockDiagonalTwoByTwoMatrixParameter
                ? new CanonicalSelectionGradientProjector.Workspace(
                        stateDimension,
                        (AbstractBlockDiagonalTwoByTwoMatrixParameter) selectionMatrixParameter)
                : null;

        branchAdjoints       = new GaussianBranchGradientAdjoints(stateDimension);
        dLogL_dFFlat         = new double[stateDimension * stateDimension];
        dLogL_dVFlat         = new double[stateDimension * stateDimension];
        rotationGradientFlat = new double[stateDimension * stateDimension];
    }

    @Override
    public boolean supportsParameter(final Parameter parameter) {
        return BlockDiagonalFormulaSupport.supportsSelectionParameter(
                selectionMatrixParameter,
                parameter);
    }

    /**
     * Computes ∂logL/∂A by iterating over branches t = 0 … T−2, computing the
     * transition-space sensitivities ∂logL/∂F_t and ∂logL/∂f_t at each branch, then
     * delegating the chain rule to {@code repr.accumulateSelectionGradient}.
     *
     * <p>The step-noise inverse Q_step_t^{-1} is computed fresh each step from a
     * Cholesky factorisation of {@code trajectory.stepCovariances[t]}.  For models with
     * a time-invariant noise covariance (e.g. uniform-grid OU process) this computation
     * could be cached; that optimisation is deferred to avoid premature complexity.
     *
     * @return flattened gradient matching the underlying matrix-parameter storage order
     */
    @Override
    public double[] computeGradient(final Parameter parameter,
                                    final BranchSmootherStats[] smootherStats,
                                    final ForwardTrajectory trajectory,
                                    final GaussianTransitionRepresentation repr,
                                    final TimeGrid timeGrid) {
        if (shouldUseBlockDiagonalNativePath(parameter, repr)) {
            return computeBlockDiagonalNativeGradient(parameter, smootherStats, trajectory, repr, timeGrid);
        }

        final int T = trajectory.timeCount;
        final int d = stateDimension;
        final double[] gradientAccumulator = new double[d * d];

        for (int t = 0; t < T - 1; ++t) {
            branchAdjoints.compute(
                    smootherStats[t],
                    smootherStats[t + 1],
                    trajectory.transitionMatrices[t],
                    trajectory.transitionOffsets[t],
                    trajectory.stepCovariances[t]);

            repr.accumulateSelectionGradient(t, t + 1, timeGrid,
                    branchAdjoints.dLogL_dF(), branchAdjoints.dLogL_df(), gradientAccumulator);

            repr.accumulateSelectionGradientFromCovariance(t, t + 1, timeGrid,
                    branchAdjoints.dLogL_dV(), gradientAccumulator);
        }

        if (parameter == selectionMatrixParameter) {
            return gradientAccumulator;
        }
        return pullBackBlockGradient(parameter, gradientAccumulator);
    }

    private boolean shouldUseBlockDiagonalNativePath(final Parameter requestedParameter,
                                                    final GaussianTransitionRepresentation repr) {
        return BlockDiagonalFormulaSupport.shouldUseNativeSelectionPath(
                selectionMatrixParameter,
                requestedParameter,
                BlockDiagonalFormulaSupport.ouProcessModel(repr));
    }

    private double[] computeBlockDiagonalNativeGradient(final Parameter requestedParameter,
                                                       final BranchSmootherStats[] smootherStats,
                                                       final ForwardTrajectory trajectory,
                                                       final GaussianTransitionRepresentation repr,
                                                       final TimeGrid timeGrid) {
        final OUProcessModel processModel = BlockDiagonalFormulaSupport.ouProcessModel(repr);
        if (processModel == null) {
            throw new IllegalArgumentException("Native block-diagonal gradient requires an OU process model");
        }
        final BlockDiagonalNativeCanonicalParameterization blockParameterization =
                BlockDiagonalFormulaSupport.nativeParameterization(processModel);
        final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter =
                (AbstractBlockDiagonalTwoByTwoMatrixParameter) selectionMatrixParameter;

        final int T = trajectory.timeCount;
        final int d = stateDimension;
        final double[] compressedDGradient = new double[blockParameter.getCompressedDDimension()];
        Arrays.fill(rotationGradientFlat, 0.0);
        final double[] stationaryMean = new double[d];
        processModel.getInitialMean(stationaryMean);

        for (int t = 0; t < T - 1; ++t) {
            final double dt = timeGrid.getDelta(t, t + 1);

            branchAdjoints.compute(
                    smootherStats[t],
                    smootherStats[t + 1],
                    trajectory.transitionMatrices[t],
                    trajectory.transitionOffsets[t],
                    trajectory.stepCovariances[t]);
            branchAdjoints.copyDLogLDFToFlat(dLogL_dFFlat);

            blockParameterization.accumulateNativeGradientFromTransitionFlat(
                    dt, stationaryMean, dLogL_dFFlat, branchAdjoints.dLogL_df(),
                    compressedDGradient, rotationGradientFlat);

            branchAdjoints.copyDLogLDVToFlat(dLogL_dVFlat);

            blockParameterization.accumulateNativeGradientFromCovarianceStationaryFlat(
                    processModel.getDiffusionMatrix(), dt, dLogL_dVFlat,
                    compressedDGradient, rotationGradientFlat);
        }

        final double[] nativeGradient = new double[blockParameter.getBlockDiagonalNParameters()];
        blockParameter.chainGradient(compressedDGradient, nativeGradient);
        return assembleBlockGradientResult(requestedParameter, blockParameter, nativeGradient, rotationGradientFlat);
    }

    private double[] assembleBlockGradientResult(final Parameter requestedParameter,
                                                 final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                                 final double[] nativeGradient,
                                                 final double[] gradientR) {
        return BlockDiagonalFormulaSupport.assembleNativeSelectionGradient(
                stateDimension,
                requestedParameter,
                blockParameter,
                nativeGradient,
                gradientR);
    }

    private double[] pullBackBlockGradient(final Parameter requestedParameter,
                                           final double[] denseGradient) {
        final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter =
                (AbstractBlockDiagonalTwoByTwoMatrixParameter) selectionMatrixParameter;
        return BlockDiagonalFormulaSupport.pullBackDenseSelectionGradient(
                stateDimension,
                requestedParameter,
                blockParameter,
                denseGradient,
                blockProjectorWorkspace);
    }

}
