package dr.evomodel.continuous.ou;

import dr.evomodel.continuous.ou.canonical.CanonicalOUKernel;
import dr.inference.model.AbstractModel;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Model;
import dr.inference.model.OrthogonalMatrixProvider;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianBranchTransitionKernel;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianBranchTransitionKernel;

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
        implements
        GaussianBranchTransitionKernel, CanonicalGaussianBranchTransitionKernel {
    /**
     * Strategy for computing the V-path contribution of ∂logL/∂A.
     */
    public enum CovarianceGradientMethod {
        VAN_LOAN_ADJOINT,
        LYAPUNOV_ADJOINT,
        STATIONARY_LYAPUNOV
    }

    private final int stateDimension;
    private final MatrixParameterInterface driftMatrix;
    private final MatrixParameterInterface diffusionMatrix;
    private final Parameter stationaryMean;
    private final MatrixParameter initialCovariance;
    private final CovarianceGradientMethod covarianceGradientMethod;
    private final SelectionMatrixParameterization selectionMatrixParameterization;
    private final CanonicalOUKernel canonicalKernel;

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
        this.stateDimension         = stateDimension;
        this.driftMatrix            = driftMatrix;
        this.diffusionMatrix        = diffusionMatrix;
        this.stationaryMean         = stationaryMean;
        this.initialCovariance      = initialCovariance;
        this.covarianceGradientMethod = covarianceGradientMethod;
        this.selectionMatrixParameterization =
                SelectionMatrixParameterizationFactory.create(driftMatrix);
        final DiffusionMatrixParameterization diffusionMatrixParameterization =
                DiffusionMatrixParameterizationFactory.create(diffusionMatrix);
        final OUCovarianceGradientStrategy covarianceGradientStrategy =
                createCovarianceGradientStrategy(covarianceGradientMethod);
        this.canonicalKernel =
                new CanonicalOUKernel(
                        stateDimension,
                        selectionMatrixParameterization,
                        diffusionMatrixParameterization,
                        diffusionMatrix,
                        covarianceGradientStrategy,
                        new CanonicalOUKernel.InitialMoments() {
                            @Override
                            public void getInitialMean(final double[] out) {
                                fillInitialMeanFromParameter(out);
                            }

                            @Override
                            public void getInitialCovarianceFlat(final double[] out) {
                                copyMatrixParameterFlat(initialCovariance, out);
                            }
                        });

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

    private static OUCovarianceGradientStrategy createCovarianceGradientStrategy(
            final CovarianceGradientMethod method) {
        switch (method) {
            case VAN_LOAN_ADJOINT:
                return new VanLoanCovarianceGradientStrategy();
            case LYAPUNOV_ADJOINT:
                return new LyapunovCovarianceGradientStrategy();
            case STATIONARY_LYAPUNOV:
                return new StationaryLyapunovCovarianceGradientStrategy();
            default:
                throw new IllegalStateException("Unknown CovarianceGradientMethod: " + method);
        }
    }

    public SelectionMatrixParameterization getSelectionMatrixParameterization() {
        return selectionMatrixParameterization;
    }

    public CanonicalOUKernel getCanonicalKernel() {
        return canonicalKernel;
    }

    @Override
    public void fillInitialCanonicalState(final CanonicalGaussianState out) {
        canonicalKernel.fillInitialCanonicalState(out);
    }

    @Override
    public void fillCanonicalTransition(final double dt, final CanonicalGaussianTransition out) {
        canonicalKernel.fillCanonicalTransition(dt, out);
    }

    @Override
    public void getInitialMean(final double[] out) {
        checkVectorLength(out, stateDimension, "initial mean");
        fillInitialMeanFromParameter(out);
    }

    @Override
    public void getInitialCovariance(final double[][] out) {
        checkSquareMatrix(out, stateDimension, "initial covariance");
        copyMatrixParameter(initialCovariance, out);
    }

    @Override
    public void fillTransitionMatrix(final double dt, final double[][] out) {
        final double[] flat = new double[stateDimension * stateDimension];
        canonicalKernel.fillTransitionMatrixFlat(dt, flat);
        copyFlatToMatrix(flat, out, stateDimension);
    }

    public void fillTransitionMatrixFlat(final double dt, final double[] out) {
        canonicalKernel.fillTransitionMatrixFlat(dt, out);
    }

    @Override
    public void fillTransitionOffset(final double dt, final double[] out) {
        canonicalKernel.fillTransitionOffset(dt, out);
    }

    @Override
    public void fillTransitionCovariance(final double dt, final double[][] out) {
        final double[] flat = new double[stateDimension * stateDimension];
        canonicalKernel.fillTransitionCovarianceFlat(dt, flat);
        copyFlatToMatrix(flat, out, stateDimension);
    }

    public void fillTransitionCovarianceFlat(final double dt, final double[] out) {
        canonicalKernel.fillTransitionCovarianceFlat(dt, out);
    }

    /**
     * Accumulates the V-path contribution of ∂logL/∂A by dispatching to the
     * strategy selected at construction time:
     * <ul>
     *   <li>{@link CovarianceGradientMethod#VAN_LOAN_ADJOINT}: backprop through the
     *       2d×2d Van Loan block-exponential (exact, one expm(2d×2d) call).</li>
     *   <li>{@link CovarianceGradientMethod#LYAPUNOV_ADJOINT}: 5-point Gauss–Legendre
     *       numerical integration of the adjoint Lyapunov ODE.</li>
     *   <li>{@link CovarianceGradientMethod#STATIONARY_LYAPUNOV}: stationary covariance
     *       Lyapunov adjoint for stable drift matrices.</li>
     * </ul>
     */
    @Override
    public void accumulateSelectionGradientFromCovariance(final double dt,
                                                          final double[][] dLogL_dV,
                                                          final double[] gradientAccumulator) {
        final double[] flat = new double[stateDimension * stateDimension];
        copyMatrixToFlat(dLogL_dV, flat, stateDimension);
        canonicalKernel.accumulateSelectionGradientFromCovarianceFlat(dt, flat, false, gradientAccumulator);
    }

    public void accumulateSelectionGradientFromCovarianceFlat(final double dt,
                                                              final double[] dLogL_dV,
                                                              final boolean transposeAdjoint,
                                                              final double[] gradientAccumulator) {
        canonicalKernel.accumulateSelectionGradientFromCovarianceFlat(
                dt, dLogL_dV, transposeAdjoint, gradientAccumulator);
    }

    @Override
    public void accumulateDiffusionGradient(final double dt,
                                            final double[][] dLogL_dV,
                                            final double[] gradientAccumulator) {
        final double[] flat = new double[stateDimension * stateDimension];
        copyMatrixToFlat(dLogL_dV, flat, stateDimension);
        canonicalKernel.accumulateDiffusionGradientFlat(dt, flat, false, gradientAccumulator);
    }

    public void accumulateDiffusionGradientFlat(final double dt,
                                                final double[] dLogL_dV,
                                                final boolean transposeAdjoint,
                                                final double[] gradientAccumulator) {
        canonicalKernel.accumulateDiffusionGradientFlat(
                dt, dLogL_dV, transposeAdjoint, gradientAccumulator);
    }

    public double contractBranchLengthGradientFlat(final double dt,
                                                   final double[] dLogL_dF,
                                                   final double[] dLogL_df,
                                                   final double[] dLogL_dV,
                                                   final boolean transposeCovarianceAdjoint) {
        return canonicalKernel.contractBranchLengthGradientFlat(
                dt,
                dLogL_dF,
                dLogL_df,
                dLogL_dV,
                transposeCovarianceAdjoint);
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
        final double[] flat = new double[stateDimension * stateDimension];
        copyMatrixToFlat(dLogL_dF, flat, stateDimension);
        canonicalKernel.accumulateSelectionGradientFlat(dt, flat, dLogL_df, gradientAccumulator);
    }

    public void accumulateSelectionGradientFlat(final double dt,
                                                final double[] dLogL_dF,
                                                final double[] dLogL_df,
                                                final double[] gradientAccumulator) {
        canonicalKernel.accumulateSelectionGradientFlat(dt, dLogL_dF, dLogL_df, gradientAccumulator);
    }

    public void accumulateSelectionGradient(final double dt,
                                            final double[] mean,
                                            final double[][] dLogL_dF,
                                            final double[] dLogL_df,
                                            final double[] gradientAccumulator) {
        final double[] flat = new double[stateDimension * stateDimension];
        copyMatrixToFlat(dLogL_dF, flat, stateDimension);
        canonicalKernel.accumulateSelectionGradientFlat(dt, mean, flat, dLogL_df, gradientAccumulator);
    }

    public void accumulateSelectionGradientFlat(final double dt,
                                                final double[] mean,
                                                final double[] dLogL_dF,
                                                final double[] dLogL_df,
                                                final double[] gradientAccumulator) {
        canonicalKernel.accumulateSelectionGradientFlat(dt, mean, dLogL_dF, dLogL_df, gradientAccumulator);
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

    private void fillInitialMeanFromParameter(final double[] out) {
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

    private static void copyMatrixParameter(final MatrixParameterInterface parameter, final double[][] out) {
        for (int i = 0; i < out.length; ++i) {
            for (int j = 0; j < out[i].length; ++j) {
                out[i][j] = parameter.getParameterValue(i, j);
            }
        }
    }

    private static void copyMatrixParameterFlat(final MatrixParameterInterface parameter, final double[] out) {
        final int dimension = parameter.getRowDimension();
        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                out[rowOffset + j] = parameter.getParameterValue(i, j);
            }
        }
    }

    private static void copyMatrixToFlat(final double[][] source, final double[] out, final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            System.arraycopy(source[i], 0, out, i * dimension, dimension);
        }
    }

    private static void copyFlatToMatrix(final double[] source, final double[][] out, final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            System.arraycopy(source, i * dimension, out[i], 0, dimension);
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
}
