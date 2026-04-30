package dr.inference.timeseries.engine.gaussian;

import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.OrthogonalBlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.OrthogonalMatrixProvider;
import dr.inference.model.Parameter;
import dr.inference.timeseries.core.TimeGrid;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockCanonicalParameterization;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;

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

    private final Parameter selectionMatrixParameter;
    private final int stateDimension;

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
        if (shouldUseOrthogonalNativePath(parameter, repr)) {
            return computeOrthogonalNativeGradient(parameter, smootherStats, trajectory, repr, timeGrid);
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

    private boolean shouldUseOrthogonalNativePath(final Parameter requestedParameter,
                                                  final GaussianTransitionRepresentation repr) {
        if (Boolean.getBoolean(DISABLE_ORTHOGONAL_NATIVE_SELECTION_PATH_PROPERTY)) {
            return false;
        }
        if (requestedParameter == selectionMatrixParameter) {
            return false;
        }
        if (!(selectionMatrixParameter instanceof OrthogonalBlockDiagonalPolarStableMatrixParameter)) {
            return false;
        }
        if (!(repr instanceof OUProcessModel)) {
            return false;
        }

        final OUProcessModel processModel = (OUProcessModel) repr;
        return processModel.getCovarianceGradientMethod() == OUProcessModel.CovarianceGradientMethod.STATIONARY_LYAPUNOV
                && processModel.getSelectionMatrixParameterization()
                instanceof OrthogonalBlockCanonicalParameterization;
    }

    private double[] computeOrthogonalNativeGradient(final Parameter requestedParameter,
                                                     final BranchSmootherStats[] smootherStats,
                                                     final ForwardTrajectory trajectory,
                                                     final GaussianTransitionRepresentation repr,
                                                     final TimeGrid timeGrid) {
        final OUProcessModel processModel = (OUProcessModel) repr;
        final OrthogonalBlockCanonicalParameterization orthogonalParameterization =
                (OrthogonalBlockCanonicalParameterization) processModel.getSelectionMatrixParameterization();
        final OrthogonalBlockDiagonalPolarStableMatrixParameter blockParameter =
                (OrthogonalBlockDiagonalPolarStableMatrixParameter) selectionMatrixParameter;

        final int T = trajectory.timeCount;
        final int d = stateDimension;
        final double[] compressedDGradient = new double[blockParameter.getCompressedDDimension()];
        final double[][] rotationGradient = new double[d][d];
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

            orthogonalParameterization.accumulateNativeGradientFromTransition(
                    dt, stationaryMean, dLogL_dF, dLogL_df,
                    compressedDGradient, rotationGradient);

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

            orthogonalParameterization.accumulateNativeGradientFromCovarianceStationary(
                    processModel.getDiffusionMatrix(), dt, dLogL_dV,
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

    private double[] pullBackBlockGradient(final Parameter requestedParameter,
                                           final double[] denseGradient) {
        final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter =
                (AbstractBlockDiagonalTwoByTwoMatrixParameter) selectionMatrixParameter;
        final int d = stateDimension;
        final double[][] gradientA = new double[d][d];
        final double[] rData = new double[d * d];
        final double[] rinvData = new double[d * d];
        final double[] rawBlockSource = new double[blockParameter.getCompressedDDimension()];
        final double[] nativeGradient = new double[blockParameter.getBlockDiagonalNParameters()];
        final double[][] gradientD = new double[d][d];
        final double[][] gradientR = new double[d][d];
        final double[][] aMatrix = new double[d][d];

        fillMatrixFromColumnMajorGradient(gradientA, denseGradient, d);
        blockParameter.fillRAndRinv(rData, rinvData);
        copyMatrixParameter(blockParameter, aMatrix);
        computeGradientWrtBlockDiagonalBasis(gradientA, rData, rinvData, gradientD);
        computeGradientWrtRotationMatrix(gradientA, aMatrix, rinvData, gradientR);
        compressActiveBlockGradient(blockParameter, gradientD, rawBlockSource);
        blockParameter.chainGradient(rawBlockSource, nativeGradient);

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

    private static void fillMatrixFromColumnMajorGradient(final double[][] out,
                                                          final double[] gradient,
                                                          final int d) {
        for (int col = 0; col < d; ++col) {
            for (int row = 0; row < d; ++row) {
                out[row][col] = gradient[col * d + row];
            }
        }
    }

    private static void computeGradientWrtBlockDiagonalBasis(final double[][] gradientA,
                                                             final double[] rData,
                                                             final double[] rinvData,
                                                             final double[][] out) {
        final int d = gradientA.length;
        final double[][] temp = new double[d][d];
        final double[][] rinvTranspose = new double[d][d];

        for (int row = 0; row < d; ++row) {
            for (int col = 0; col < d; ++col) {
                rinvTranspose[row][col] = rinvData[col * d + row];
            }
        }

        multiplyTransposeLeft(rData, d, gradientA, temp);
        multiply(temp, rinvTranspose, out);
    }

    private static void computeGradientWrtRotationMatrix(final double[][] gradientA,
                                                         final double[][] aMatrix,
                                                         final double[] rinvData,
                                                         final double[][] out) {
        final int d = gradientA.length;
        final double[][] aTranspose = new double[d][d];
        final double[][] temp = new double[d][d];
        final double[][] rinvTranspose = new double[d][d];

        transpose(aMatrix, aTranspose);
        multiply(gradientA, aTranspose, temp);
        multiply(aTranspose, gradientA, out);

        for (int row = 0; row < d; ++row) {
            for (int col = 0; col < d; ++col) {
                rinvTranspose[row][col] = rinvData[col * d + row];
            }
        }

        subtractInPlace(temp, out);
        multiply(temp, rinvTranspose, out);
    }

    private static void compressActiveBlockGradient(final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                                    final double[][] gradientD,
                                                    final double[] out) {
        final int d = gradientD.length;
        int blockIndex = 0;
        final int upperBase = d;
        final int lowerBase = d + blockParameter.getNum2x2Blocks();

        for (int i = 0; i < d; ++i) {
            out[i] = gradientD[i][i];
        }

        for (int b = 0; b < blockParameter.getNumBlocks(); ++b) {
            if (blockParameter.getBlockSizes()[b] != 2) {
                continue;
            }
            final int start = blockParameter.getBlockStarts()[b];
            out[upperBase + blockIndex] = gradientD[start][start + 1];
            out[lowerBase + blockIndex] = gradientD[start + 1][start];
            blockIndex++;
        }
    }

    private static void copyMatrixParameter(final MatrixParameter parameter, final double[][] out) {
        final int d = out.length;
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                out[i][j] = parameter.getParameterValue(i, j);
            }
        }
    }

    private static void transpose(final double[][] in, final double[][] out) {
        final int d = in.length;
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                out[i][j] = in[j][i];
            }
        }
    }

    private static void multiplyTransposeLeft(final double[] leftRowMajor,
                                              final int d,
                                              final double[][] right,
                                              final double[][] out) {
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += leftRowMajor[k * d + i] * right[k][j];
                }
                out[i][j] = sum;
            }
        }
    }

    private static void multiply(final double[][] left,
                                 final double[][] right,
                                 final double[][] out) {
        final int d = left.length;
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += left[i][k] * right[k][j];
                }
                out[i][j] = sum;
            }
        }
    }

    private static void subtractInPlace(final double[][] minuend,
                                        final double[][] subtrahend) {
        final int d = minuend.length;
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                minuend[i][j] -= subtrahend[i][j];
            }
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
