package dr.evomodel.continuous.ou;

import dr.inference.model.AbstractModel;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Model;
import dr.inference.model.OrthogonalMatrixProvider;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.timeseries.core.LatentProcessModel;
import dr.inference.timeseries.core.TimeGrid;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianBranchTransitionKernel;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianUtils;
import dr.evomodel.treedatalikelihood.continuous.gaussian.GaussianBranchTransitionKernel;
import dr.inference.timeseries.gaussian.DiffusionMatrixParameterization;
import dr.inference.timeseries.gaussian.DiffusionMatrixParameterizationFactory;
import dr.evomodel.continuous.ou.SelectionMatrixParameterization;
import dr.evomodel.continuous.ou.SelectionMatrixParameterizationFactory;
import dr.inference.timeseries.representation.GaussianComputationMode;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;
import dr.inference.timeseries.representation.KernelBackedGaussianTransitionRepresentation;
import dr.inference.timeseries.representation.RepresentableProcess;

/**
 * Multivariate Ornstein-Uhlenbeck process model with exact Gaussian transitions.
 *
 * <p>The latent diffusion is
 * <pre>
 *   dX_t = -A (X_t - mu) dt + dW_t,
 *   E[dW_t dW_t^T] = Q dt.
 * </pre>
 * For a time increment dt, the exact discrete-time transition is
 * <pre>
 *   X_{t+dt} | X_t ~ N(F(dt) X_t + b(dt), V(dt)),
 *   F(dt) = exp(-A dt),
 *   b(dt) = (I - F(dt)) mu,
 *   V(dt) = integral_0^dt exp(-A s) Q exp(-A^T s) ds.
 * </pre>
 *
 * <p>The transition covariance is computed by the standard Van Loan block-exponential
 * construction.  The selection-matrix gradient hook uses <em>adjoint/backpropagation</em>
 * through the exact transition map rather than looping over forward basis directions.
 *
 * <h3>V-path gradient strategies</h3>
 * The gradient ∂logL/∂A has two contributions: the <em>F/f-path</em> (through the
 * transition matrix and offset) and the <em>V-path</em> (through the step covariance).
 * Two strategies for the V-path are available via {@link CovarianceGradientMethod}:
 * <ul>
 *   <li>{@link CovarianceGradientMethod#VAN_LOAN_ADJOINT}: backpropagates
 *       through the 2d×2d Van Loan block-exponential using the Fréchet-adjoint identity.
 *       Exact, O(d³) per branch after a single expm(2d×2d) call.</li>
 *   <li>{@link CovarianceGradientMethod#LYAPUNOV_ADJOINT}: numerically integrates the
 *       adjoint of the Lyapunov ODE via 5-point Gauss–Legendre quadrature:
 *       <pre>∂logL/∂A ≈ -2 ∫₀^{dt} Ψ(s) V(s) ds</pre>
 *       where Ψ(s) = F(dt−s)ᵀ G_V F(dt−s) and V(s) = Van Loan at step s.
 *       Mathematically equivalent to VAN_LOAN_ADJOINT; useful as an independent check.</li>
 *   <li>{@link CovarianceGradientMethod#STATIONARY_LYAPUNOV}: uses the stationary-covariance
 *       identity <pre>V(dt) = S - F(dt) S F(dt)^T</pre> where {@code S} solves
 *       <pre>A S + S A^T = Q</pre>. The covariance-path gradient is split into the
 *       stationary-covariance adjoint and the F-inside-V adjoint. This requires a stable
 *       drift matrix; zero-drift and other non-stable cases are rejected.</li>
 * </ul>
 * By default, orthogonal block-diagonal drift parametrizations use
 * {@link CovarianceGradientMethod#STATIONARY_LYAPUNOV}, while general dense drift matrices
 * use {@link CovarianceGradientMethod#VAN_LOAN_ADJOINT}. The active strategy is selected
 * at construction time and is immutable unless explicitly overridden.
 */
public class OUProcessModel extends AbstractModel
        implements LatentProcessModel, RepresentableProcess,
        GaussianBranchTransitionKernel, CanonicalGaussianBranchTransitionKernel {
    private static final String DIFFUSION_GRADIENT_QUADRATURE_SUBSTEPS_PROPERTY =
            "beast.experimental.ouDiffusionGradientQuadratureSubsteps";

    /**
     * Strategy for computing the V-path contribution of ∂logL/∂A.
     */
//    public enum CovarianceGradientMethod {
//        /**
//         * Backpropagates through the 2d×2d Van Loan block-exponential using the
//         * Fréchet-adjoint identity. Exact and efficient (one expm(2d×2d) call per branch).
//         */
//        VAN_LOAN_ADJOINT,
//        /**
//         * Numerically integrates the adjoint Lyapunov ODE via 5-point Gauss–Legendre
//         * quadrature: ∂logL/∂A ≈ -2∫₀^{dt} Ψ(s)V(s) ds, where
//         * Ψ(s) = F(dt−s)ᵀ G_V F(dt−s).  Mathematically equivalent to VAN_LOAN_ADJOINT;
//         * more expensive but useful as an independent numerical check.
//         */
//        LYAPUNOV_ADJOINT
//    }
    public enum CovarianceGradientMethod {
        VAN_LOAN_ADJOINT,
        LYAPUNOV_ADJOINT,
        STATIONARY_LYAPUNOV
    }

    // ── 5-point Gauss–Legendre nodes and weights on [-1, 1] ─────────────────────────
    private static final double[] GL5_NODES = {
            0.0,
            -0.5384693101056831,  0.5384693101056831,
            -0.9061798459386640,  0.9061798459386640
    };
    private static final double[] GL5_WEIGHTS = {
            0.5688888888888889,
            0.4786286704993665,  0.4786286704993665,
            0.2369268850561891,  0.2369268850561891
    };

    private final int stateDimension;
    private final MatrixParameterInterface driftMatrix;
    private final MatrixParameterInterface diffusionMatrix;
    private final Parameter stationaryMean;
    private final MatrixParameter initialCovariance;
    private final CovarianceGradientMethod covarianceGradientMethod;
    private final GaussianComputationMode defaultComputationMode;
    private final SelectionMatrixParameterization selectionMatrixParameterization;
    private final DiffusionMatrixParameterization diffusionMatrixParameterization;
    private final GaussianTransitionRepresentation transitionRepresentation;
    private final ThreadLocal<Workspace> workspaceLocal;

    /**
     * Constructs an OUProcessModel using the default covariance-adjoint strategy for the
     * supplied drift parametrization:
     * orthogonal block-diagonal drift matrices use stationary-Lyapunov adjoints, while
     * general dense drift matrices use the Van Loan adjoint.
     */
    public OUProcessModel(final String name,
                          final int stateDimension,
                          final MatrixParameterInterface driftMatrix,
                          final MatrixParameterInterface diffusionMatrix,
                          final Parameter stationaryMean,
                          final MatrixParameter initialCovariance) {
        this(name, stateDimension, driftMatrix, diffusionMatrix, stationaryMean,
                initialCovariance, defaultCovarianceGradientMethod(driftMatrix));
    }

    /** Constructs an OUProcessModel with an explicit V-path gradient strategy. */
    public OUProcessModel(final String name,
                          final int stateDimension,
                          final MatrixParameterInterface driftMatrix,
                          final MatrixParameterInterface diffusionMatrix,
                          final Parameter stationaryMean,
                          final MatrixParameter initialCovariance,
                          final CovarianceGradientMethod covarianceGradientMethod) {
        this(name, stateDimension, driftMatrix, diffusionMatrix, stationaryMean,
                initialCovariance, covarianceGradientMethod, GaussianComputationMode.EXPECTATION);
    }

    public OUProcessModel(final String name,
                          final int stateDimension,
                          final MatrixParameterInterface driftMatrix,
                          final MatrixParameterInterface diffusionMatrix,
                          final Parameter stationaryMean,
                          final MatrixParameter initialCovariance,
                          final CovarianceGradientMethod covarianceGradientMethod,
                          final GaussianComputationMode defaultComputationMode) {
        super(name);
        if (stateDimension < 1) {
            throw new IllegalArgumentException("stateDimension must be at least 1");
        }
        if (driftMatrix == null) {
            throw new IllegalArgumentException("driftMatrix must not be null");
        }
        if (diffusionMatrix == null) {
            throw new IllegalArgumentException("diffusionMatrix must not be null");
        }
        if (stationaryMean == null) {
            throw new IllegalArgumentException("stationaryMean must not be null");
        }
        if (initialCovariance == null) {
            throw new IllegalArgumentException("initialCovariance must not be null");
        }
        if (covarianceGradientMethod == null) {
            throw new IllegalArgumentException("covarianceGradientMethod must not be null");
        }
        if (defaultComputationMode == null) {
            throw new IllegalArgumentException("defaultComputationMode must not be null");
        }
        this.stateDimension         = stateDimension;
        this.driftMatrix            = driftMatrix;
        this.diffusionMatrix        = diffusionMatrix;
        this.stationaryMean         = stationaryMean;
        this.initialCovariance      = initialCovariance;
        this.covarianceGradientMethod = covarianceGradientMethod;
        this.defaultComputationMode = defaultComputationMode;
        this.selectionMatrixParameterization =
                SelectionMatrixParameterizationFactory.create(driftMatrix);
        this.diffusionMatrixParameterization =
                DiffusionMatrixParameterizationFactory.create(diffusionMatrix);
        this.transitionRepresentation =
                new KernelBackedGaussianTransitionRepresentation(this);
        final boolean allocateStationaryLyapunov =
                covarianceGradientMethod == CovarianceGradientMethod.STATIONARY_LYAPUNOV;
        this.workspaceLocal = new ThreadLocal<Workspace>() {
            @Override
            protected Workspace initialValue() {
                return new Workspace(stateDimension, allocateStationaryLyapunov);
            }
        };

        validateShapes();

        addVariable(driftMatrix);
        addVariable(diffusionMatrix);
        addVariable(stationaryMean);
        addVariable(initialCovariance);
    }

    @Override
    public int getStateDimension() {
        return stateDimension;
    }

    public MatrixParameterInterface getDriftMatrix() {
        return driftMatrix;
    }

    public MatrixParameterInterface getDiffusionMatrix() {
        return diffusionMatrix;
    }

    public Parameter getStationaryMeanParameter() {
        return stationaryMean;
    }

    public MatrixParameter getInitialCovarianceParameter() {
        return initialCovariance;
    }

    public CovarianceGradientMethod getCovarianceGradientMethod() {
        return covarianceGradientMethod;
    }

    public static boolean usesOrthogonalBlockSelectionChart(
            final MatrixParameterInterface driftMatrix) {
        if (driftMatrix instanceof AbstractBlockDiagonalTwoByTwoMatrixParameter) {
            final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter =
                    (AbstractBlockDiagonalTwoByTwoMatrixParameter) driftMatrix;
            return blockParameter.getRotationMatrixParameter() instanceof OrthogonalMatrixProvider;
        }
        return false;
    }

    public static CovarianceGradientMethod defaultCovarianceGradientMethod(
            final MatrixParameterInterface driftMatrix) {
        if (usesOrthogonalBlockSelectionChart(driftMatrix)) {
            return CovarianceGradientMethod.STATIONARY_LYAPUNOV;
        }
        return CovarianceGradientMethod.VAN_LOAN_ADJOINT;
    }

    public GaussianComputationMode getDefaultComputationMode() {
        return defaultComputationMode;
    }

    public SelectionMatrixParameterization getSelectionMatrixParameterization() {
        return selectionMatrixParameterization;
    }

    @Override
    public <T> boolean supportsRepresentation(final Class<T> representationClass) {
        return representationClass.isAssignableFrom(GaussianTransitionRepresentation.class)
                || representationClass.isAssignableFrom(GaussianBranchTransitionKernel.class)
                || representationClass.isAssignableFrom(CanonicalGaussianBranchTransitionKernel.class);
    }

    @Override
    @SuppressWarnings("unchecked")
    public <T> T getRepresentation(final Class<T> representationClass) {
        if (!supportsRepresentation(representationClass)) {
            throw new IllegalArgumentException("Unsupported representation: " + representationClass.getName());
        }
        if (representationClass.isAssignableFrom(GaussianTransitionRepresentation.class)) {
            return (T) transitionRepresentation;
        }
        return (T) this;
    }

    /**
     * Backward-compatible transition accessor that delegates to the current
     * representation API.
     */
    public void getTransitionMatrix(final int fromIndex,
                                    final int toIndex,
                                    final TimeGrid timeGrid,
                                    final double[][] out) {
        transitionRepresentation.getTransitionMatrix(fromIndex, toIndex, timeGrid, out);
    }

    /**
     * Backward-compatible transition accessor that delegates to the current
     * representation API.
     */
    public void getTransitionOffset(final int fromIndex,
                                    final int toIndex,
                                    final TimeGrid timeGrid,
                                    final double[] out) {
        transitionRepresentation.getTransitionOffset(fromIndex, toIndex, timeGrid, out);
    }

    /**
     * Backward-compatible transition accessor that delegates to the current
     * representation API.
     */
    public void getTransitionCovariance(final int fromIndex,
                                        final int toIndex,
                                        final TimeGrid timeGrid,
                                        final double[][] out) {
        transitionRepresentation.getTransitionCovariance(fromIndex, toIndex, timeGrid, out);
    }

    /**
     * Backward-compatible index-based gradient hook routed through the current
     * representation API.
     */
    public void accumulateSelectionGradient(final int fromIndex,
                                            final int toIndex,
                                            final TimeGrid timeGrid,
                                            final double[][] dLogL_dF,
                                            final double[] dLogL_df,
                                            final double[] gradientAccumulator) {
        transitionRepresentation.accumulateSelectionGradient(
                fromIndex, toIndex, timeGrid, dLogL_dF, dLogL_df, gradientAccumulator);
    }

    /**
     * Backward-compatible index-based covariance gradient hook routed through
     * the current representation API.
     */
    public void accumulateSelectionGradientFromCovariance(final int fromIndex,
                                                          final int toIndex,
                                                          final TimeGrid timeGrid,
                                                          final double[][] dLogL_dV,
                                                          final double[] gradientAccumulator) {
        transitionRepresentation.accumulateSelectionGradientFromCovariance(
                fromIndex, toIndex, timeGrid, dLogL_dV, gradientAccumulator);
    }

    @Override
    public void fillInitialCanonicalState(final CanonicalGaussianState out) {
        final Workspace workspace = workspace();
        final double[] mean = workspace.vector0;
        final double[][] covariance = workspace.squareMatrices[0];
        getInitialMean(mean);
        getInitialCovariance(covariance);
        CanonicalGaussianUtils.fillStateFromMoments(mean, covariance, out);
    }

    @Override
    public void fillCanonicalTransition(final double dt, final CanonicalGaussianTransition out) {
        final Workspace workspace = workspace();
        if (selectionMatrixParameterization instanceof OrthogonalBlockDiagonalSelectionMatrixParameterization) {
            final double[] mean = workspace.vector0;
            getInitialMean(mean);
            ((OrthogonalBlockDiagonalSelectionMatrixParameterization) selectionMatrixParameterization)
                    .fillCanonicalTransition(diffusionMatrix, mean, dt, out);
            return;
        }
        final double[][] transitionMatrix = workspace.squareMatrices[0];
        final double[] transitionOffset = workspace.vector1;
        final double[][] transitionCovariance = workspace.squareMatrices[1];
        fillTransitionMatrix(dt, transitionMatrix);
        fillTransitionOffset(dt, transitionOffset);
        fillTransitionCovariance(dt, transitionCovariance);
        CanonicalGaussianUtils.fillTransitionFromMoments(
                transitionMatrix,
                transitionOffset,
                transitionCovariance,
                out);
    }

    @Override
    public void getInitialMean(final double[] out) {
        checkVectorLength(out, stateDimension, "initial mean");
        if (stationaryMean.getDimension() == 1) {
            final double value = stationaryMean.getParameterValue(0);
            for (int i = 0; i < out.length; ++i) {
                out[i] = value;
            }
        } else if (stationaryMean.getDimension() == stateDimension) {
            for (int i = 0; i < out.length; ++i) {
                out[i] = stationaryMean.getParameterValue(i);
            }
        } else {
            throw new IllegalStateException("stationaryMean dimension must be 1 or equal to the state dimension");
        }
    }

    @Override
    public void getInitialCovariance(final double[][] out) {
        checkSquareMatrix(out, stateDimension, "initial covariance");
        copyMatrixParameter(initialCovariance, out);
    }

    @Override
    public void fillTransitionMatrix(final double dt, final double[][] out) {
        checkSquareMatrix(out, stateDimension, "transition matrix");
        selectionMatrixParameterization.fillTransitionMatrix(dt, out);
    }

    public void fillTransitionMatrixFlat(final double dt, final double[] out) {
        checkFlatSquare(out, stateDimension, "transition matrix");
        selectionMatrixParameterization.fillTransitionMatrixFlat(dt, out);
    }

    @Override
    public void fillTransitionOffset(final double dt, final double[] out) {
        checkVectorLength(out, stateDimension, "transition offset");
        final Workspace workspace = workspace();
        final double[][] transitionMatrix = workspace.squareMatrices[2];
        final double[] mu = workspace.vector0;
        final double[] transformedMean = workspace.vector2;

        fillTransitionMatrix(dt, transitionMatrix);
        getInitialMean(mu);
        MatrixExponentialUtils.multiply(transitionMatrix, mu, transformedMean);

        for (int i = 0; i < stateDimension; ++i) {
            out[i] = mu[i] - transformedMean[i];
        }
    }

    @Override
    public void fillTransitionCovariance(final double dt, final double[][] out) {
        checkSquareMatrix(out, stateDimension, "transition covariance");
        final Workspace workspace = workspace();

        if (selectionMatrixParameterization instanceof OrthogonalBlockDiagonalSelectionMatrixParameterization) {
            ((OrthogonalBlockDiagonalSelectionMatrixParameterization) selectionMatrixParameterization)
                    .fillTransitionCovariance(diffusionMatrix, dt, out);
            return;
        }

        final double[][] a = workspace.squareMatrices[2];
        final double[][] q = workspace.squareMatrices[3];
        selectionMatrixParameterization.fillSelectionMatrix(a);
        diffusionMatrixParameterization.fillDiffusionMatrix(q);

        final int blockDimension = 2 * stateDimension;
        final double[][] vanLoan = workspace.blockMatrices[0];
        final double[][] vanLoanExp = workspace.blockMatrices[1];

        for (int i = 0; i < stateDimension; ++i) {
            for (int j = 0; j < stateDimension; ++j) {
                vanLoan[i][j] = -dt * a[i][j];
                vanLoan[i][j + stateDimension] = dt * q[i][j];
                vanLoan[i + stateDimension][j] = 0.0;
                vanLoan[i + stateDimension][j + stateDimension] = dt * a[j][i];
            }
        }

        MatrixExponentialUtils.expm(vanLoan, vanLoanExp);

        for (int i = 0; i < stateDimension; ++i) {
            for (int j = 0; j < stateDimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < stateDimension; ++k) {
                    sum += vanLoanExp[i][k + stateDimension] * vanLoanExp[j][k];
                }
                out[i][j] = sum;
            }
        }
        MatrixExponentialUtils.symmetrize(out);
    }

    /**
     * Accumulates the V-path contribution of ∂logL/∂A by dispatching to the
     * strategy selected at construction time:
     * <ul>
     *   <li>{@link CovarianceGradientMethod#VAN_LOAN_ADJOINT}: backprop through the
     *       2d×2d Van Loan block-exponential (exact, one expm(2d×2d) call).</li>
     *   <li>{@link CovarianceGradientMethod#LYAPUNOV_ADJOINT}: 5-point Gauss–Legendre
     *       numerical integration of the adjoint Lyapunov ODE.</li>
     * </ul>
     */
    @Override
    public void accumulateSelectionGradientFromCovariance(final double dt,
                                                          final double[][] dLogL_dV,
                                                          final double[] gradientAccumulator) {
        final int d = stateDimension;
        final Workspace workspace = workspace();

        // Read current parameter values into local arrays (shared by both strategies)
        final double[][] a = workspace.squareMatrices[0];
        final double[][] q = workspace.squareMatrices[1];
        selectionMatrixParameterization.fillSelectionMatrix(a);
        diffusionMatrixParameterization.fillDiffusionMatrix(q);

        switch (covarianceGradientMethod) {
            case VAN_LOAN_ADJOINT:
                accumulateViaVanLoanAdjoint(dt, d, a, q, dLogL_dV, gradientAccumulator);
                break;
            case LYAPUNOV_ADJOINT:
                accumulateViaLyapunovAdjoint(dt, d, a, q, dLogL_dV, gradientAccumulator);
                break;
            case STATIONARY_LYAPUNOV:
                accumulateViaStationaryLyapunov(dt, d, a, q, dLogL_dV, gradientAccumulator);
                break;
            default:
                throw new IllegalStateException("Unknown CovarianceGradientMethod: " + covarianceGradientMethod);
        }
    }

    public void accumulateSelectionGradientFromCovarianceFlat(final double dt,
                                                              final double[] dLogL_dV,
                                                              final boolean transposeAdjoint,
                                                              final double[] gradientAccumulator) {
        checkFlatSquare(dLogL_dV, stateDimension, "covariance adjoint");
        final int d = stateDimension;
        final Workspace workspace = workspace();

        final double[][] a = workspace.squareMatrices[0];
        final double[][] q = workspace.squareMatrices[1];
        selectionMatrixParameterization.fillSelectionMatrix(a);
        diffusionMatrixParameterization.fillDiffusionMatrix(q);

        switch (covarianceGradientMethod) {
            case VAN_LOAN_ADJOINT:
                accumulateViaVanLoanAdjointFlat(dt, d, a, q, dLogL_dV, transposeAdjoint, gradientAccumulator);
                break;
            case LYAPUNOV_ADJOINT:
                accumulateViaLyapunovAdjointFlat(dt, d, a, q, dLogL_dV, transposeAdjoint, gradientAccumulator);
                break;
            case STATIONARY_LYAPUNOV:
                accumulateViaStationaryLyapunovFlat(dt, d, a, q, dLogL_dV, transposeAdjoint, gradientAccumulator);
                break;
            default:
                throw new IllegalStateException("Unknown CovarianceGradientMethod: " + covarianceGradientMethod);
        }
    }

    @Override
    public void accumulateDiffusionGradient(final double dt,
                                            final double[][] dLogL_dV,
                                            final double[] gradientAccumulator) {
        final int d = stateDimension;
        final Workspace workspace = workspace();
        final double[][] a = workspace.squareMatrices[0];
        final double[][] gSym = workspace.squareMatrices[2];
        final double[][] fS = workspace.squareMatrices[3];
        final double[][] fST = workspace.squareMatrices[4];
        final double[][] tempDxD = workspace.squareMatrices[5];
        final double[][] contrib = workspace.squareMatrices[6];
        final double[][] expScratch = workspace.squareMatrices[7];

        selectionMatrixParameterization.fillSelectionMatrix(a);
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                gSym[i][j] = 0.5 * (dLogL_dV[i][j] + dLogL_dV[j][i]);
            }
        }

        final int substeps = Math.max(1,
                Integer.getInteger(DIFFUSION_GRADIENT_QUADRATURE_SUBSTEPS_PROPERTY, 2));
        final double h = dt / substeps;
        for (int sub = 0; sub < substeps; ++sub) {
            final double t0 = sub * h;
            for (int idx = 0; idx < GL5_NODES.length; ++idx) {
                final double s = t0 + 0.5 * h * (GL5_NODES[idx] + 1.0);
                final double scaledWeight = 0.5 * h * GL5_WEIGHTS[idx];

                buildExpmMinusAs(s, a, d, fS, expScratch);
                MatrixExponentialUtils.transpose(fS, fST);
                MatrixExponentialUtils.multiply(fST, gSym, tempDxD);
                MatrixExponentialUtils.multiply(tempDxD, fS, contrib);

                for (int i = 0; i < d; ++i) {
                    for (int j = 0; j < d; ++j) {
                        gradientAccumulator[i * d + j] += scaledWeight * contrib[i][j];
                    }
                }
            }
        }
    }

    public void accumulateDiffusionGradientFlat(final double dt,
                                                final double[] dLogL_dV,
                                                final boolean transposeAdjoint,
                                                final double[] gradientAccumulator) {
        checkFlatSquare(dLogL_dV, stateDimension, "covariance adjoint");
        final int d = stateDimension;
        final Workspace workspace = workspace();
        final double[][] a = workspace.squareMatrices[0];
        final double[][] gSym = workspace.squareMatrices[2];
        final double[][] fS = workspace.squareMatrices[3];
        final double[][] fST = workspace.squareMatrices[4];
        final double[][] tempDxD = workspace.squareMatrices[5];
        final double[][] contrib = workspace.squareMatrices[6];
        final double[][] expScratch = workspace.squareMatrices[7];

        selectionMatrixParameterization.fillSelectionMatrix(a);
        fillSymmetricFromFlat(dLogL_dV, transposeAdjoint, d, gSym);

        final int substeps = Math.max(1,
                Integer.getInteger(DIFFUSION_GRADIENT_QUADRATURE_SUBSTEPS_PROPERTY, 2));
        final double h = dt / substeps;
        for (int sub = 0; sub < substeps; ++sub) {
            final double t0 = sub * h;
            for (int idx = 0; idx < GL5_NODES.length; ++idx) {
                final double s = t0 + 0.5 * h * (GL5_NODES[idx] + 1.0);
                final double scaledWeight = 0.5 * h * GL5_WEIGHTS[idx];

                buildExpmMinusAs(s, a, d, fS, expScratch);
                MatrixExponentialUtils.transpose(fS, fST);
                MatrixExponentialUtils.multiply(fST, gSym, tempDxD);
                MatrixExponentialUtils.multiply(tempDxD, fS, contrib);

                for (int i = 0; i < d; ++i) {
                    for (int j = 0; j < d; ++j) {
                        gradientAccumulator[i * d + j] += scaledWeight * contrib[i][j];
                    }
                }
            }
        }
    }

    /**
     * V-path gradient via Van Loan adjoint backpropagation.
     *
     * <p>The Van Loan block and its exponential are:
     * <pre>
     *   M   = [[-A dt,  Q dt], [0,  A^T dt]]
     *   E   = expm(M) = [[F,  V·F^{-T}], [0,  F^{-T}]]
     * </pre>
     * so  V = E_{01} · F^T.  Differentiating:
     * <pre>
     *   G_E = [[  0,       G_V · F  ],
     *          [  0,   −V · G_V · F ]]
     *   G_M = adjointExp(M, G_E)
     *   ∂logL/∂A_{kl} += −dt · G_M[k][l]  +  dt · G_M[d+l][d+k]
     * </pre>
     */
    private void accumulateViaVanLoanAdjoint(final double dt, final int d,
                                              final double[][] a, final double[][] q,
                                              final double[][] dLogL_dV,
                                              final double[] gradientAccumulator) {
        final Workspace workspace = workspace();
        final int blockDim = 2 * d;

        // Build Van Loan block M = [[-A dt, Q dt], [0, A^T dt]]
        final double[][] vanLoanMatrix = workspace.blockMatrices[0];
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                vanLoanMatrix[i][j]         = -dt * a[i][j];
                vanLoanMatrix[i][j + d]     =  dt * q[i][j];
                vanLoanMatrix[i + d][j]     =  0.0;
                vanLoanMatrix[i + d][j + d] =  dt * a[j][i];   // A^T_{ij} = A_{ji}
            }
        }

        final double[][] vanLoanExp = workspace.blockMatrices[1];
        MatrixExponentialUtils.expm(vanLoanMatrix, vanLoanExp);

        // GV_F = G_V · F = G_V · E_{00}
        final double[][] gvF = workspace.squareMatrices[2];
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += dLogL_dV[i][k] * vanLoanExp[k][j];
                }
                gvF[i][j] = sum;
            }
        }

        // Recover V = E_{01} · F^T  (= E_{01} · E_{00}^T)
        final double[][] v = workspace.squareMatrices[3];
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += vanLoanExp[i][k + d] * vanLoanExp[j][k];
                }
                v[i][j] = sum;
            }
        }

        // Upstream G_E on the 2d×2d block
        final double[][] upstreamBlock = workspace.blockMatrices[2];
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                upstreamBlock[i][j + d] = gvF[i][j];             // top-right: G_V · F
                double vgvf = 0.0;
                for (int k = 0; k < d; ++k) {
                    vgvf += v[i][k] * gvF[k][j];
                }
                upstreamBlock[i + d][j + d] = -vgvf;             // bot-right: −V · G_V · F
            }
        }

        final double[][] gradM = workspace.blockMatrices[3];
        MatrixExponentialUtils.adjointExp(vanLoanMatrix, upstreamBlock, gradM);

        for (int k = 0; k < d; ++k) {
            for (int l = 0; l < d; ++l) {
                gradientAccumulator[k * d + l] +=
                        -dt * gradM[k][l] + dt * gradM[d + l][d + k];
            }
        }
    }

    private void accumulateViaVanLoanAdjointFlat(final double dt, final int d,
                                                 final double[][] a, final double[][] q,
                                                 final double[] dLogL_dV,
                                                 final boolean transposeAdjoint,
                                                 final double[] gradientAccumulator) {
        final Workspace workspace = workspace();

        final double[][] vanLoanMatrix = workspace.blockMatrices[0];
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                vanLoanMatrix[i][j] = -dt * a[i][j];
                vanLoanMatrix[i][j + d] = dt * q[i][j];
                vanLoanMatrix[i + d][j] = 0.0;
                vanLoanMatrix[i + d][j + d] = dt * a[j][i];
            }
        }

        final double[][] vanLoanExp = workspace.blockMatrices[1];
        MatrixExponentialUtils.expm(vanLoanMatrix, vanLoanExp);

        final double[][] gvF = workspace.squareMatrices[2];
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += flatSquareValue(dLogL_dV, i, k, d, transposeAdjoint) * vanLoanExp[k][j];
                }
                gvF[i][j] = sum;
            }
        }

        final double[][] v = workspace.squareMatrices[3];
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += vanLoanExp[i][k + d] * vanLoanExp[j][k];
                }
                v[i][j] = sum;
            }
        }

        final double[][] upstreamBlock = workspace.blockMatrices[2];
        clearSquare(upstreamBlock);
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                upstreamBlock[i][j + d] = gvF[i][j];
                double vgvf = 0.0;
                for (int k = 0; k < d; ++k) {
                    vgvf += v[i][k] * gvF[k][j];
                }
                upstreamBlock[i + d][j + d] = -vgvf;
            }
        }

        final double[][] gradM = workspace.blockMatrices[3];
        MatrixExponentialUtils.adjointExp(vanLoanMatrix, upstreamBlock, gradM);

        for (int k = 0; k < d; ++k) {
            for (int l = 0; l < d; ++l) {
                gradientAccumulator[k * d + l] +=
                        -dt * gradM[k][l] + dt * gradM[d + l][d + k];
            }
        }
    }

    private void accumulateViaLyapunovAdjoint(final double dt, final int d,
                                              final double[][] a, final double[][] q,
                                              final double[][] dLogL_dV,
                                              final double[] gradientAccumulator) {
        final Workspace workspace = workspace();
        final double[][] gSym = workspace.squareMatrices[2];
        final double[][] fRemaining = workspace.squareMatrices[3];
        final double[][] fRemainingT = workspace.squareMatrices[4];
        final double[][] psi = workspace.squareMatrices[5];
        final double[][] tempDxD = workspace.squareMatrices[6];
        final double[][] vS = workspace.squareMatrices[7];
        final double[][] expScratch = workspace.squareMatrices[8];

        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                gSym[i][j] = 0.5 * (dLogL_dV[i][j] + dLogL_dV[j][i]);
            }
        }

        for (int idx = 0; idx < GL5_NODES.length; ++idx) {
            final double s = 0.5 * dt * (GL5_NODES[idx] + 1.0);
            final double scaledWeight = 0.5 * dt * GL5_WEIGHTS[idx];

            buildExpmMinusAs(dt - s, a, d, fRemaining, expScratch);
            MatrixExponentialUtils.transpose(fRemaining, fRemainingT);
            MatrixExponentialUtils.multiply(fRemainingT, gSym, tempDxD);
            MatrixExponentialUtils.multiply(tempDxD, fRemaining, psi);
            buildVanLoanCovariance(s, a, q, d, vS, workspace.blockMatrices[0], workspace.blockMatrices[1]);

            MatrixExponentialUtils.multiply(vS, psi, tempDxD);
            for (int i = 0; i < d; ++i) {
                for (int j = 0; j < d; ++j) {
                    gradientAccumulator[i * d + j] += -2.0 * scaledWeight * tempDxD[i][j];
                }
            }
        }
    }

    private void accumulateViaLyapunovAdjointFlat(final double dt, final int d,
                                                  final double[][] a, final double[][] q,
                                                  final double[] dLogL_dV,
                                                  final boolean transposeAdjoint,
                                                  final double[] gradientAccumulator) {
        final Workspace workspace = workspace();
        final double[][] gSym = workspace.squareMatrices[2];
        final double[][] fRemaining = workspace.squareMatrices[3];
        final double[][] fRemainingT = workspace.squareMatrices[4];
        final double[][] psi = workspace.squareMatrices[5];
        final double[][] tempDxD = workspace.squareMatrices[6];
        final double[][] vS = workspace.squareMatrices[7];
        final double[][] expScratch = workspace.squareMatrices[8];

        fillSymmetricFromFlat(dLogL_dV, transposeAdjoint, d, gSym);

        for (int idx = 0; idx < GL5_NODES.length; ++idx) {
            final double s = 0.5 * dt * (GL5_NODES[idx] + 1.0);
            final double scaledWeight = 0.5 * dt * GL5_WEIGHTS[idx];

            buildExpmMinusAs(dt - s, a, d, fRemaining, expScratch);
            MatrixExponentialUtils.transpose(fRemaining, fRemainingT);
            MatrixExponentialUtils.multiply(fRemainingT, gSym, tempDxD);
            MatrixExponentialUtils.multiply(tempDxD, fRemaining, psi);
            buildVanLoanCovariance(s, a, q, d, vS, workspace.blockMatrices[0], workspace.blockMatrices[1]);

            MatrixExponentialUtils.multiply(vS, psi, tempDxD);
            for (int i = 0; i < d; ++i) {
                for (int j = 0; j < d; ++j) {
                    gradientAccumulator[i * d + j] += -2.0 * scaledWeight * tempDxD[i][j];
                }
            }
        }
    }

    private void accumulateViaStationaryLyapunov(final double dt, final int d,
                                                 final double[][] a, final double[][] q,
                                                 final double[][] dLogL_dV,
                                                 final double[] gradientAccumulator) {
        final Workspace workspace = workspace();
        final double[][] gV = workspace.squareMatrices[4];

        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                gV[i][j] = 0.5 * (dLogL_dV[i][j] + dLogL_dV[j][i]);
            }
        }
        accumulateViaStationaryLyapunovSymmetric(dt, d, a, q, gV, gradientAccumulator);
    }

    private void accumulateViaStationaryLyapunovFlat(final double dt, final int d,
                                                     final double[][] a, final double[][] q,
                                                     final double[] dLogL_dV,
                                                     final boolean transposeAdjoint,
                                                     final double[] gradientAccumulator) {
        final Workspace workspace = workspace();
        final double[][] gV = workspace.squareMatrices[4];
        fillSymmetricFromFlat(dLogL_dV, transposeAdjoint, d, gV);
        accumulateViaStationaryLyapunovSymmetric(dt, d, a, q, gV, gradientAccumulator);
    }

    private void accumulateViaStationaryLyapunovSymmetric(final double dt, final int d,
                                                          final double[][] a, final double[][] q,
                                                          final double[][] gV,
                                                          final double[] gradientAccumulator) {
        final Workspace workspace = workspace();
        final double[][] sStat = workspace.squareMatrices[2];
        final double[][] f = workspace.squareMatrices[3];
        final double[][] tempDxD = workspace.squareMatrices[5];
        final double[][] gS = workspace.squareMatrices[6];
        final double[][] y = workspace.squareMatrices[7];
        final double[][] yT = workspace.squareMatrices[8];
        final double[][] gAStationary = workspace.squareMatrices[9];
        final double[][] sStatT = workspace.squareMatrices[10];
        final double[][] gFCov = workspace.squareMatrices[11];
        final double[][] gVT = workspace.squareMatrices[12];
        final double[][] tempDxD2 = workspace.squareMatrices[13];
        final double[][] gX = workspace.squareMatrices[14];
        final double[][] fT = workspace.squareMatrices[15];
        final double[][] expScratch = workspace.squareMatrices[16];
        final double[][] minusAdt = workspace.squareMatrices[17];

        solveStationaryLyapunovDense(a, q, d, sStat, workspace);
        buildExpmMinusAs(dt, a, d, f, expScratch);

        // G_S = G_V - F^T G_V F
        MatrixExponentialUtils.transpose(f, fT);
        MatrixExponentialUtils.multiply(fT, gV, tempDxD);
        MatrixExponentialUtils.multiply(tempDxD, f, gS);
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                gS[i][j] = gV[i][j] - gS[i][j];
            }
        }

        solveAdjointStationaryLyapunovDense(a, gS, d, y, workspace);

        // G_A^(S) = -(Y S + Y^T S)
        MatrixExponentialUtils.transpose(y, yT);
        MatrixExponentialUtils.multiply(y, sStat, gAStationary);
        MatrixExponentialUtils.multiply(yT, sStat, tempDxD);
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                gAStationary[i][j] = -(gAStationary[i][j] + tempDxD[i][j]);
            }
        }

        // G_F^(cov) = -(G_V F S + G_V^T F S^T)
        MatrixExponentialUtils.transpose(sStat, sStatT);
        MatrixExponentialUtils.multiply(gV, f, tempDxD);
        MatrixExponentialUtils.multiply(tempDxD, sStat, gFCov);
        MatrixExponentialUtils.transpose(gV, gVT);
        MatrixExponentialUtils.multiply(gVT, f, tempDxD2);
        MatrixExponentialUtils.multiply(tempDxD2, sStatT, tempDxD);
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                gFCov[i][j] = -(gFCov[i][j] + tempDxD[i][j]);
            }
        }

        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                minusAdt[i][j] = -dt * a[i][j];
            }
        }
        MatrixExponentialUtils.adjointExp(minusAdt, gFCov, gX);

        // Shared accumulator uses row-major flattening of the raw [d×d] contribution.
        // The stationary covariance path contributes the transpose of
        // (G_A^(S) - dt * G_X^T) to match the parameter-indexing convention used
        // by the rest of the analytical engine.
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                gradientAccumulator[i * d + j] += gAStationary[j][i] - dt * gX[i][j];
            }
        }
    }

    private static void solveStationaryLyapunovDense(final double[][] a,
                                                     final double[][] q,
                                                     final int d,
                                                     final double[][] sOut,
                                                     final Workspace workspace) {
        final int n = d * d;
        final double[][] augmented = workspace.lyapunovAugmented;
        if (augmented == null) {
            throw new IllegalStateException("Stationary Lyapunov workspace is not available");
        }

        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                final int row = i * d + j;
                for (int p = 0; p < d; ++p) {
                    for (int r = 0; r < d; ++r) {
                        final int col = p * d + r;
                        double value = 0.0;
                        if (i == p) {
                            value += a[j][r];
                        }
                        if (j == r) {
                            value += a[i][p];
                        }
                        augmented[row][col] = value;
                    }
                }
                augmented[row][n] = q[i][j];
            }
        }

        try {
            solveAugmentedLinearSystemInPlace(augmented);
        } catch (IllegalStateException e) {
            throw new IllegalArgumentException(
                    "Stationary covariance is only defined for stable drift matrices", e);
        }

        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                sOut[i][j] = augmented[i * d + j][n];
            }
        }
        MatrixExponentialUtils.symmetrize(sOut);
    }

    private static void solveAdjointStationaryLyapunovDense(final double[][] a,
                                                            final double[][] gS,
                                                            final int d,
                                                            final double[][] yOut,
                                                            final Workspace workspace) {
        final double[][] aT = workspace.squareMatrices[17];
        MatrixExponentialUtils.transpose(a, aT);
        solveStationaryLyapunovDense(aT, gS, d, yOut, workspace);
    }

    private static void solveAugmentedLinearSystemInPlace(final double[][] augmented) {
        final int n = augmented.length;
        for (int col = 0; col < n; ++col) {
            int pivot = col;
            double maxAbs = Math.abs(augmented[col][col]);
            for (int row = col + 1; row < n; ++row) {
                final double abs = Math.abs(augmented[row][col]);
                if (abs > maxAbs) {
                    maxAbs = abs;
                    pivot = row;
                }
            }

            if (maxAbs <= 1e-14) {
                throw new IllegalStateException("Singular system");
            }

            if (pivot != col) {
                final double[] tmp = augmented[col];
                augmented[col] = augmented[pivot];
                augmented[pivot] = tmp;
            }

            final double diag = augmented[col][col];
            for (int j = col; j <= n; ++j) {
                augmented[col][j] /= diag;
            }

            for (int row = 0; row < n; ++row) {
                if (row == col) continue;
                final double factor = augmented[row][col];
                if (factor == 0.0) continue;
                for (int j = col; j <= n; ++j) {
                    augmented[row][j] -= factor * augmented[col][j];
                }
            }
        }
    }

    /** Computes out = expm(-A * t). */
    private static void buildExpmMinusAs(final double t, final double[][] a,
                                         final int d,
                                         final double[][] out,
                                         final double[][] scaledMinusA) {
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                scaledMinusA[i][j] = -t * a[i][j];
            }
        }
        MatrixExponentialUtils.expm(scaledMinusA, out);
    }

    /**
     * Computes out = V(t) = ∫₀ᵗ exp(−Au) Q exp(−Aᵀu) du via the Van Loan construction.
     */
    private static void buildVanLoanCovariance(final double t, final double[][] a,
                                               final double[][] q, final int d,
                                               final double[][] out,
                                               final double[][] vanLoan,
                                               final double[][] exp) {
        final int blockDim = 2 * d;
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                vanLoan[i][j]         = -t * a[i][j];
                vanLoan[i][j + d]     =  t * q[i][j];
                vanLoan[i + d][j]     =  0.0;
                vanLoan[i + d][j + d] =  t * a[j][i];  // A^T
            }
        }
        MatrixExponentialUtils.expm(vanLoan, exp);
        // V = E_{01} · F^T = E_{01} · E_{00}^T
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += exp[i][k + d] * exp[j][k];
                }
                out[i][j] = sum;
            }
        }
    }

    /**
     * Backpropagates through the exact transition map
     * <pre>
     *   X = -dt A,
     *   F = exp(X),
     *   b = (I - F) mu.
     * </pre>
     *
     * <p>Let G_F = dL/dF and g_b = dL/db be the transition-space sensitivities provided
     * by the smoother-based gradient formula.  Since
     * <pre>
     *   db = -(dF) mu,
     * </pre>
     * the total upstream sensitivity on F is
     * <pre>
     *   G_total = G_F - g_b mu^T.
     * </pre>
     * We then backpropagate through F = exp(X) using the adjoint Fr\'echet identity
     * and finally through X = -dt A.
     *
     * <p>No forward basis-direction loop is used here: the implementation computes a
     * single matrix adjoint and accumulates it directly into the flattened A-gradient.
     */
    @Override
    public void accumulateSelectionGradient(final double dt,
                                            final double[][] dLogL_dF,
                                            final double[] dLogL_df,
                                            final double[] gradientAccumulator) {
        final Workspace workspace = workspace();
        final double[] mu = workspace.vector0;
        getInitialMean(mu);
        accumulateSelectionGradient(dt, mu, dLogL_dF, dLogL_df, gradientAccumulator);
    }

    public void accumulateSelectionGradientFlat(final double dt,
                                                final double[] dLogL_dF,
                                                final double[] dLogL_df,
                                                final double[] gradientAccumulator) {
        checkFlatSquare(dLogL_dF, stateDimension, "transition adjoint");
        final Workspace workspace = workspace();
        final double[] mu = workspace.vector0;
        getInitialMean(mu);
        accumulateSelectionGradientFlat(dt, mu, dLogL_dF, dLogL_df, gradientAccumulator);
    }

    public void accumulateSelectionGradient(final double dt,
                                            final double[] mean,
                                            final double[][] dLogL_dF,
                                            final double[] dLogL_df,
                                            final double[] gradientAccumulator) {
        selectionMatrixParameterization.accumulateGradientFromTransition(
                dt, mean, dLogL_dF, dLogL_df, gradientAccumulator);
    }

    public void accumulateSelectionGradientFlat(final double dt,
                                                final double[] mean,
                                                final double[] dLogL_dF,
                                                final double[] dLogL_df,
                                                final double[] gradientAccumulator) {
        checkFlatSquare(dLogL_dF, stateDimension, "transition adjoint");
        selectionMatrixParameterization.accumulateGradientFromTransitionFlat(
                dt, mean, dLogL_dF, dLogL_df, gradientAccumulator);
    }

    @Override
    protected void handleModelChangedEvent(final Model model, final Object object, final int index) {
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(final Variable variable,
                                              final int index,
                                              final Parameter.ChangeType type) {
        fireModelChanged();
    }

    @Override
    protected void storeState() {
        // no-op
    }

    @Override
    protected void restoreState() {
        // no-op
    }

    @Override
    protected void acceptState() {
        // no-op
    }

    private Workspace workspace() {
        return workspaceLocal.get();
    }

    private void validateShapes() {
        if (driftMatrix.getRowDimension() != stateDimension || driftMatrix.getColumnDimension() != stateDimension) {
            throw new IllegalArgumentException("driftMatrix must be stateDimension x stateDimension");
        }
        if (diffusionMatrix.getRowDimension() != stateDimension || diffusionMatrix.getColumnDimension() != stateDimension) {
            throw new IllegalArgumentException("diffusionMatrix must be stateDimension x stateDimension");
        }
        if (initialCovariance.getRowDimension() != stateDimension || initialCovariance.getColumnDimension() != stateDimension) {
            throw new IllegalArgumentException("initialCovariance must be stateDimension x stateDimension");
        }
        if (stationaryMean.getDimension() != 1 && stationaryMean.getDimension() != stateDimension) {
            throw new IllegalArgumentException("stationaryMean dimension must be 1 or stateDimension");
        }
    }

    private static double validatedDelta(final TimeGrid timeGrid, final int fromIndex, final int toIndex) {
        final double dt = timeGrid.getDelta(fromIndex, toIndex);
        if (!(dt > 0.0)) {
            throw new IllegalArgumentException("Time increments must be strictly positive");
        }
        return dt;
    }

    private static void copyMatrixParameter(final MatrixParameterInterface parameter, final double[][] out) {
        for (int i = 0; i < out.length; ++i) {
            for (int j = 0; j < out[i].length; ++j) {
                out[i][j] = parameter.getParameterValue(i, j);
            }
        }
    }

    private static void checkVectorLength(final double[] vector, final int expectedLength, final String label) {
        if (vector == null || vector.length != expectedLength) {
            throw new IllegalArgumentException(label + " must have length " + expectedLength);
        }
    }

    private static void checkSquareMatrix(final double[][] matrix, final int expectedSize, final String label) {
        if (matrix == null || matrix.length != expectedSize) {
            throw new IllegalArgumentException(label + " must be a " + expectedSize + "x" + expectedSize + " matrix");
        }
        for (double[] row : matrix) {
            if (row == null || row.length != expectedSize) {
                throw new IllegalArgumentException(label + " must be a " + expectedSize + "x" + expectedSize + " matrix");
            }
        }
    }

    private static void checkFlatSquare(final double[] matrix, final int expectedSize, final String label) {
        if (matrix == null || matrix.length != expectedSize * expectedSize) {
            throw new IllegalArgumentException(
                    label + " must have length " + (expectedSize * expectedSize)
                            + " for a " + expectedSize + "x" + expectedSize + " matrix");
        }
    }

    private static double flatSquareValue(final double[] source,
                                          final int row,
                                          final int col,
                                          final int dim,
                                          final boolean transpose) {
        return transpose ? source[col * dim + row] : source[row * dim + col];
    }

    private static void fillSymmetricFromFlat(final double[] source,
                                              final boolean transpose,
                                              final int dim,
                                              final double[][] out) {
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                out[i][j] = 0.5 * (
                        flatSquareValue(source, i, j, dim, transpose)
                                + flatSquareValue(source, j, i, dim, transpose));
            }
        }
    }

    private static void clearSquare(final double[][] matrix) {
        for (double[] row : matrix) {
            for (int j = 0; j < row.length; ++j) {
                row[j] = 0.0;
            }
        }
    }

    private static final class Workspace {
        private static final int SQUARE_MATRIX_COUNT = 18;
        private static final int BLOCK_MATRIX_COUNT = 4;

        private final double[][][] squareMatrices;
        private final double[][][] blockMatrices;
        private final double[][] flatAdapterMatrix;
        private final double[] vector0;
        private final double[] vector1;
        private final double[] vector2;
        private final double[][] lyapunovAugmented;

        private Workspace(final int dimension,
                          final boolean allocateStationaryLyapunov) {
            this.squareMatrices = allocateMatrixStack(SQUARE_MATRIX_COUNT, dimension, dimension);
            this.blockMatrices = allocateMatrixStack(BLOCK_MATRIX_COUNT, 2 * dimension, 2 * dimension);
            this.flatAdapterMatrix = allocateMatrix(dimension, dimension);
            this.vector0 = new double[dimension];
            this.vector1 = new double[dimension];
            this.vector2 = new double[dimension];
            this.lyapunovAugmented = allocateStationaryLyapunov
                    ? allocateMatrix(dimension * dimension, dimension * dimension + 1)
                    : null;
        }
    }

    private static double[][][] allocateMatrixStack(final int count,
                                                    final int rows,
                                                    final int cols) {
        final double[][][] stack = new double[count][][];
        for (int i = 0; i < count; ++i) {
            stack[i] = allocateMatrix(rows, cols);
        }
        return stack;
    }

    private static double[][] allocateMatrix(final int rows,
                                             final int cols) {
        final double[][] matrix = new double[rows][];
        for (int i = 0; i < rows; ++i) {
            matrix[i] = new double[cols];
        }
        return matrix;
    }
}
