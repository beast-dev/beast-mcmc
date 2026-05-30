package dr.inference.timeseries.engine.gaussian;

import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.OrthogonalMatrixProvider;
import dr.inference.model.Parameter;
import dr.evomodel.treedatalikelihood.continuous.canonical.gradient.CanonicalSelectionGradientProjector;
import dr.inference.timeseries.core.TimeGrid;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.BlockDiagonalNativeCanonicalParameterization;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;
import dr.inference.timeseries.representation.KernelBackedGaussianTransitionRepresentation;

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

    private static final String DISABLE_ORTHOGONAL_NATIVE_SELECTION_PATH_PROPERTY =
            "beast.experimental.disableOrthogonalNativeSelectionPath";
    private static final String DISABLE_BLOCK_DIAGONAL_NATIVE_SELECTION_PATH_PROPERTY =
            "beast.experimental.disableBlockDiagonalNativeSelectionPath";

    private final Parameter selectionMatrixParameter;
    private final int stateDimension;
    private final CanonicalSelectionGradientProjector.Workspace blockProjectorWorkspace;

    // Pre-allocated working arrays — reused across computeGradient calls.
    private final double[]   meanResidual;          // r̄_t
    private final double[][] crossCov;              // V_{t+1,t} = P_{t+1|T} · G_t^T
    private final double[][] branchMat;             // B_t = V_{t+1,t} − F_t P_{t|T} + r̄_t m_t^T
    private final double[][] stepCovInv;            // V_t^{-1}  (copied per step, then inverted)
    private final double[][] tempDxD;               // scratch [d × d]
    private final double[][] tempDxD2;              // scratch [d × d]
    private final double[]   tempD;                 // scratch [d]
    private final double[][] dLogL_dF;              // ∂logL/∂F_t
    private final double[]   dLogL_df;              // ∂logL/∂f_t
    private final double[][] residualSecondMoment;  // S_t = E[δ_t δ_t^T | Y]
    private final double[][] dLogL_dV;              // ∂logL/∂V_t = ½(V^{-1} S V^{-1} − V^{-1})
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

        meanResidual         = new double[stateDimension];
        crossCov             = new double[stateDimension][stateDimension];
        branchMat            = new double[stateDimension][stateDimension];
        stepCovInv           = new double[stateDimension][stateDimension];
        tempDxD              = new double[stateDimension][stateDimension];
        tempDxD2             = new double[stateDimension][stateDimension];
        tempD                = new double[stateDimension];
        dLogL_dF             = new double[stateDimension][stateDimension];
        dLogL_df             = new double[stateDimension];
        residualSecondMoment = new double[stateDimension][stateDimension];
        dLogL_dV             = new double[stateDimension][stateDimension];
        dLogL_dFFlat         = new double[stateDimension * stateDimension];
        dLogL_dVFlat         = new double[stateDimension * stateDimension];
        rotationGradientFlat = new double[stateDimension * stateDimension];
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
            final BranchSmootherStats curr = smootherStats[t];
            final BranchSmootherStats next = smootherStats[t + 1];
            final double[][] F_t = trajectory.transitionMatrices[t];
            final double[] f_t   = trajectory.transitionOffsets[t];

            // ── V_{t+1,t} = P_{t+1|T} · G_t^T ──────────────────────────────────
            KalmanLikelihoodEngine.multiplyMatrixMatrixTransposedRight(
                    next.smoothedCovariance, curr.smootherGain, crossCov);

            // ── r̄_t = m_{t+1|T} − F_t · m_{t|T} − f_t ─────────────────────────
            KalmanLikelihoodEngine.multiplyMatrixVector(F_t, curr.smoothedMean, meanResidual);
            for (int i = 0; i < d; ++i) {
                meanResidual[i] = next.smoothedMean[i] - meanResidual[i] - f_t[i];
            }

            // ── B_t = V_{t+1,t} − F_t · P_{t|T} + r̄_t · m_{t|T}^T ─────────────
            // Start from cross-covariance
            KalmanLikelihoodEngine.copyMatrix(crossCov, branchMat);

            // Subtract F_t · P_{t|T}
            KalmanLikelihoodEngine.multiplyMatrixMatrix(F_t, curr.smoothedCovariance, tempDxD);
            for (int i = 0; i < d; ++i) {
                for (int j = 0; j < d; ++j) {
                    branchMat[i][j] -= tempDxD[i][j];
                }
            }

            // Add outer product r̄_t · m_{t|T}^T
            for (int i = 0; i < d; ++i) {
                for (int j = 0; j < d; ++j) {
                    branchMat[i][j] += meanResidual[i] * curr.smoothedMean[j];
                }
            }

            // ── Q_step_t^{-1} ────────────────────────────────────────────────────
            KalmanLikelihoodEngine.copyMatrix(trajectory.stepCovariances[t], stepCovInv);
            final KalmanLikelihoodEngine.CholeskyFactor stepChol =
                    KalmanLikelihoodEngine.cholesky(stepCovInv);
            KalmanLikelihoodEngine.invertPositiveDefiniteFromCholesky(stepCovInv, stepChol);

            // ── ∂logL/∂F_t = Q_step_t^{-1} · B_t ───────────────────────────────
            KalmanLikelihoodEngine.multiplyMatrixMatrix(stepCovInv, branchMat, dLogL_dF);

            // ── ∂logL/∂f_t = Q_step_t^{-1} · r̄_t ──────────────────────────────
            KalmanLikelihoodEngine.multiplyMatrixVector(stepCovInv, meanResidual, dLogL_df);

            // ── Chain rule (F/f-path) → accumulate ∂logL/∂A ────────────────────
            repr.accumulateSelectionGradient(t, t + 1, timeGrid,
                    dLogL_dF, dLogL_df, gradientAccumulator);

            // ── V-path: ∂logL/∂V_t = ½(V_t^{-1} S_t V_t^{-1} − V_t^{-1}) ─────
            //
            // S_t = E[(x_{t+1} − F_t x_t − f_t)(x_{t+1} − F_t x_t − f_t)^T | Y]
            //     = P_{t+1|T}  −  F_t V_{t+1,t}^T  −  V_{t+1,t} F_t^T
            //                  +  F_t P_{t|T} F_t^T  +  r̄_t r̄_t^T

            // Start from P_{t+1|T}
            KalmanLikelihoodEngine.copyMatrix(next.smoothedCovariance, residualSecondMoment);

            // Subtract F_t × V_{t+1,t}^T  (= F_t × crossCov^T)
            KalmanLikelihoodEngine.multiplyMatrixMatrixTransposedRight(F_t, crossCov, tempDxD);
            for (int i = 0; i < d; ++i) {
                for (int j = 0; j < d; ++j) {
                    residualSecondMoment[i][j] -= tempDxD[i][j];
                }
            }

            // Subtract V_{t+1,t} × F_t^T  (= crossCov × F_t^T)
            KalmanLikelihoodEngine.multiplyMatrixMatrixTransposedRight(crossCov, F_t, tempDxD);
            for (int i = 0; i < d; ++i) {
                for (int j = 0; j < d; ++j) {
                    residualSecondMoment[i][j] -= tempDxD[i][j];
                }
            }

            // Add F_t P_{t|T} F_t^T
            KalmanLikelihoodEngine.multiplyMatrixMatrix(F_t, curr.smoothedCovariance, tempDxD);
            KalmanLikelihoodEngine.multiplyMatrixMatrixTransposedRight(tempDxD, F_t, tempDxD2);
            for (int i = 0; i < d; ++i) {
                for (int j = 0; j < d; ++j) {
                    residualSecondMoment[i][j] += tempDxD2[i][j];
                }
            }

            // Add r̄_t r̄_t^T
            for (int i = 0; i < d; ++i) {
                for (int j = 0; j < d; ++j) {
                    residualSecondMoment[i][j] += meanResidual[i] * meanResidual[j];
                }
            }

            // dLogL_dV = ½ (V_t^{-1} S_t V_t^{-1} − V_t^{-1})
            //          = ½ V_t^{-1} (S_t − V_t) V_t^{-1}
            // stepCovInv already holds V_t^{-1} from the F/f-path computation above.
            KalmanLikelihoodEngine.multiplyMatrixMatrix(stepCovInv, residualSecondMoment, tempDxD);
            KalmanLikelihoodEngine.multiplyMatrixMatrix(tempDxD, stepCovInv, dLogL_dV);
            for (int i = 0; i < d; ++i) {
                for (int j = 0; j < d; ++j) {
                    dLogL_dV[i][j] = 0.5 * (dLogL_dV[i][j] - stepCovInv[i][j]);
                }
            }

            // ── Chain rule (V-path) → accumulate ∂logL/∂A ──────────────────────
            repr.accumulateSelectionGradientFromCovariance(t, t + 1, timeGrid,
                    dLogL_dV, gradientAccumulator);
        }

        if (parameter == selectionMatrixParameter) {
            return gradientAccumulator;
        }
        return pullBackBlockGradient(parameter, gradientAccumulator);
    }

    private boolean shouldUseBlockDiagonalNativePath(final Parameter requestedParameter,
                                                    final GaussianTransitionRepresentation repr) {
        if (Boolean.getBoolean(DISABLE_BLOCK_DIAGONAL_NATIVE_SELECTION_PATH_PROPERTY)
                || Boolean.getBoolean(DISABLE_ORTHOGONAL_NATIVE_SELECTION_PATH_PROPERTY)) {
            return false;
        }
        if (requestedParameter == selectionMatrixParameter) {
            return false;
        }
        if (!(selectionMatrixParameter instanceof AbstractBlockDiagonalTwoByTwoMatrixParameter)) {
            return false;
        }
        final OUProcessModel processModel = ouProcessModel(repr);
        if (processModel == null) {
            return false;
        }
        return processModel.getCovarianceGradientMethod() == OUProcessModel.CovarianceGradientMethod.STATIONARY_LYAPUNOV
                && processModel.getSelectionMatrixParameterization()
                instanceof BlockDiagonalNativeCanonicalParameterization;
    }

    private double[] computeBlockDiagonalNativeGradient(final Parameter requestedParameter,
                                                       final BranchSmootherStats[] smootherStats,
                                                       final ForwardTrajectory trajectory,
                                                       final GaussianTransitionRepresentation repr,
                                                       final TimeGrid timeGrid) {
        final OUProcessModel processModel = ouProcessModel(repr);
        if (processModel == null) {
            throw new IllegalArgumentException("Native block-diagonal gradient requires an OU process model");
        }
        final BlockDiagonalNativeCanonicalParameterization blockParameterization =
                (BlockDiagonalNativeCanonicalParameterization) processModel.getSelectionMatrixParameterization();
        final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter =
                (AbstractBlockDiagonalTwoByTwoMatrixParameter) selectionMatrixParameter;

        final int T = trajectory.timeCount;
        final int d = stateDimension;
        final double[] compressedDGradient = new double[blockParameter.getCompressedDDimension()];
        Arrays.fill(rotationGradientFlat, 0.0);
        final double[] stationaryMean = new double[d];
        processModel.getInitialMean(stationaryMean);

        for (int t = 0; t < T - 1; ++t) {
            final BranchSmootherStats curr = smootherStats[t];
            final BranchSmootherStats next = smootherStats[t + 1];
            final double[][] F_t = trajectory.transitionMatrices[t];
            final double[] f_t   = trajectory.transitionOffsets[t];
            final double dt = timeGrid.getDelta(t, t + 1);

            KalmanLikelihoodEngine.multiplyMatrixMatrixTransposedRight(
                    next.smoothedCovariance, curr.smootherGain, crossCov);

            KalmanLikelihoodEngine.multiplyMatrixVector(F_t, curr.smoothedMean, meanResidual);
            for (int i = 0; i < d; ++i) {
                meanResidual[i] = next.smoothedMean[i] - meanResidual[i] - f_t[i];
            }

            KalmanLikelihoodEngine.copyMatrix(crossCov, branchMat);
            KalmanLikelihoodEngine.multiplyMatrixMatrix(F_t, curr.smoothedCovariance, tempDxD);
            for (int i = 0; i < d; ++i) {
                for (int j = 0; j < d; ++j) {
                    branchMat[i][j] -= tempDxD[i][j];
                    branchMat[i][j] += meanResidual[i] * curr.smoothedMean[j];
                }
            }

            KalmanLikelihoodEngine.copyMatrix(trajectory.stepCovariances[t], stepCovInv);
            final KalmanLikelihoodEngine.CholeskyFactor stepChol =
                    KalmanLikelihoodEngine.cholesky(stepCovInv);
            KalmanLikelihoodEngine.invertPositiveDefiniteFromCholesky(stepCovInv, stepChol);

            KalmanLikelihoodEngine.multiplyMatrixMatrix(stepCovInv, branchMat, dLogL_dF);
            KalmanLikelihoodEngine.multiplyMatrixVector(stepCovInv, meanResidual, dLogL_df);
            copyMatrixToFlat(dLogL_dF, dLogL_dFFlat, d);

            blockParameterization.accumulateNativeGradientFromTransitionFlat(
                    dt, stationaryMean, dLogL_dFFlat, dLogL_df,
                    compressedDGradient, rotationGradientFlat);

            KalmanLikelihoodEngine.copyMatrix(next.smoothedCovariance, residualSecondMoment);
            KalmanLikelihoodEngine.multiplyMatrixMatrixTransposedRight(F_t, crossCov, tempDxD);
            for (int i = 0; i < d; ++i) {
                for (int j = 0; j < d; ++j) {
                    residualSecondMoment[i][j] -= tempDxD[i][j];
                }
            }
            KalmanLikelihoodEngine.multiplyMatrixMatrixTransposedRight(crossCov, F_t, tempDxD);
            for (int i = 0; i < d; ++i) {
                for (int j = 0; j < d; ++j) {
                    residualSecondMoment[i][j] -= tempDxD[i][j];
                }
            }
            KalmanLikelihoodEngine.multiplyMatrixMatrix(F_t, curr.smoothedCovariance, tempDxD);
            KalmanLikelihoodEngine.multiplyMatrixMatrixTransposedRight(tempDxD, F_t, tempDxD2);
            for (int i = 0; i < d; ++i) {
                for (int j = 0; j < d; ++j) {
                    residualSecondMoment[i][j] += tempDxD2[i][j];
                    residualSecondMoment[i][j] += meanResidual[i] * meanResidual[j];
                }
            }

            KalmanLikelihoodEngine.multiplyMatrixMatrix(stepCovInv, residualSecondMoment, tempDxD);
            KalmanLikelihoodEngine.multiplyMatrixMatrix(tempDxD, stepCovInv, dLogL_dV);
            for (int i = 0; i < d; ++i) {
                for (int j = 0; j < d; ++j) {
                    dLogL_dV[i][j] = 0.5 * (dLogL_dV[i][j] - stepCovInv[i][j]);
                }
            }
            copyMatrixToFlat(dLogL_dV, dLogL_dVFlat, d);

            blockParameterization.accumulateNativeGradientFromCovarianceStationaryFlat(
                    processModel.getDiffusionMatrix(), dt, dLogL_dVFlat,
                    compressedDGradient, rotationGradientFlat);
        }

        final double[] nativeGradient = new double[blockParameter.getBlockDiagonalNParameters()];
        blockParameter.chainGradient(compressedDGradient, nativeGradient);
        return assembleBlockGradientResult(requestedParameter, blockParameter, nativeGradient, rotationGradientFlat);
    }

    private static OUProcessModel ouProcessModel(final GaussianTransitionRepresentation repr) {
        if (repr instanceof OUProcessModel) {
            return (OUProcessModel) repr;
        }
        if (repr instanceof KernelBackedGaussianTransitionRepresentation) {
            return ((KernelBackedGaussianTransitionRepresentation) repr).getProcessModel();
        }
        return null;
    }

    private double[] assembleBlockGradientResult(final Parameter requestedParameter,
                                                 final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                                 final double[] nativeGradient,
                                                 final double[] gradientR) {
        return CanonicalSelectionGradientProjector.assembleBlockGradientResultFlat(
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
        return CanonicalSelectionGradientProjector.pullBackDenseGradientToBlock(
                stateDimension,
                requestedParameter,
                blockParameter,
                denseGradient,
                blockProjectorWorkspace);
    }

    private static void copyMatrixToFlat(final double[][] source,
                                         final double[] out,
                                         final int d) {
        for (int row = 0; row < d; ++row) {
            System.arraycopy(source[row], 0, out, row * d, d);
        }
    }

}
