package dr.inference.timeseries.gaussian;

import dr.inference.model.AbstractModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.timeseries.core.LatentProcessModel;
import dr.inference.timeseries.core.TimeGrid;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianBranchTransitionKernel;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianUtils;
import dr.evomodel.treedatalikelihood.continuous.gaussian.GaussianBranchTransitionKernel;
import dr.inference.timeseries.representation.GaussianComputationMode;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;
import dr.inference.timeseries.representation.KernelBackedGaussianTransitionRepresentation;
import dr.inference.timeseries.representation.RepresentableProcess;

/**
 * First-order Euler approximation of the multivariate Ornstein-Uhlenbeck process.
 *
 * <p>The transition parameters are approximated as:
 * <pre>
 *   F = I − dt·A
 *   b = dt·A·μ
 *   V = dt·Q
 * </pre>
 * These are O(dt) accurate.  For production use prefer {@link OUProcessModel}, which
 * provides the exact matrix-exponential transitions.  This class is retained as a
 * reference implementation and for gradient cross-checks at very small dt.
 */
public class EulerOUProcessModel extends AbstractModel
        implements LatentProcessModel, RepresentableProcess,
        GaussianBranchTransitionKernel, CanonicalGaussianBranchTransitionKernel {

    private final int stateDimension;
    private final MatrixParameter driftMatrix;
    private final MatrixParameterInterface diffusionMatrix;
    private final Parameter stationaryMean;
    private final MatrixParameter initialCovariance;
    private final GaussianComputationMode defaultComputationMode;
    private final DiffusionMatrixParameterization diffusionMatrixParameterization;
    private final GaussianTransitionRepresentation transitionRepresentation;

    public EulerOUProcessModel(final String name,
                               final int stateDimension,
                               final MatrixParameter driftMatrix,
                               final MatrixParameterInterface diffusionMatrix,
                               final Parameter stationaryMean,
                               final MatrixParameter initialCovariance) {
        this(name, stateDimension, driftMatrix, diffusionMatrix,
                stationaryMean, initialCovariance, GaussianComputationMode.EXPECTATION);
    }

    public EulerOUProcessModel(final String name,
                               final int stateDimension,
                               final MatrixParameter driftMatrix,
                               final MatrixParameterInterface diffusionMatrix,
                               final Parameter stationaryMean,
                               final MatrixParameter initialCovariance,
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
        if (defaultComputationMode == null) {
            throw new IllegalArgumentException("defaultComputationMode must not be null");
        }
        this.stateDimension   = stateDimension;
        this.driftMatrix      = driftMatrix;
        this.diffusionMatrix  = diffusionMatrix;
        this.stationaryMean   = stationaryMean;
        this.initialCovariance = initialCovariance;
        this.defaultComputationMode = defaultComputationMode;
        this.diffusionMatrixParameterization =
                DiffusionMatrixParameterizationFactory.create(diffusionMatrix);
        this.transitionRepresentation =
                new KernelBackedGaussianTransitionRepresentation(this);

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

    public MatrixParameter getDriftMatrix() {
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

    public GaussianComputationMode getDefaultComputationMode() {
        return defaultComputationMode;
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

    @Override
    public void fillInitialCanonicalState(final CanonicalGaussianState out) {
        final double[] mean = new double[stateDimension];
        final double[][] covariance = new double[stateDimension][stateDimension];
        getInitialMean(mean);
        getInitialCovariance(covariance);
        CanonicalGaussianUtils.fillStateFromMoments(mean, covariance, out);
    }

    @Override
    public void fillCanonicalTransition(final double dt, final CanonicalGaussianTransition out) {
        final double[][] transitionMatrix = new double[stateDimension][stateDimension];
        final double[] transitionOffset = new double[stateDimension];
        final double[][] transitionCovariance = new double[stateDimension][stateDimension];
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
            throw new IllegalStateException(
                    "stationaryMean dimension must be 1 or equal to the state dimension");
        }
    }

    @Override
    public void getInitialCovariance(final double[][] out) {
        checkSquareMatrix(out, stateDimension, "initial covariance");
        copyMatrixParameter(initialCovariance, out);
    }

    /** F = I − dt·A */
    @Override
    public void fillTransitionMatrix(final double dt, final double[][] out) {
        checkSquareMatrix(out, stateDimension, "transition matrix");
        setIdentity(out);
        for (int i = 0; i < stateDimension; ++i) {
            for (int j = 0; j < stateDimension; ++j) {
                out[i][j] -= dt * driftMatrix.getParameterValue(i, j);
            }
        }
    }

    /** b = dt·A·μ */
    @Override
    public void fillTransitionOffset(final double dt, final double[] out) {
        checkVectorLength(out, stateDimension, "transition offset");
        final double[] mu = new double[stateDimension];
        getInitialMean(mu);
        for (int i = 0; i < stateDimension; ++i) {
            double driftTimesMean = 0.0;
            for (int j = 0; j < stateDimension; ++j) {
                driftTimesMean += driftMatrix.getParameterValue(i, j) * mu[j];
            }
            out[i] = dt * driftTimesMean;
        }
    }

    /** V = dt·Q */
    @Override
    public void fillTransitionCovariance(final double dt, final double[][] out) {
        checkSquareMatrix(out, stateDimension, "transition covariance");
        diffusionMatrixParameterization.fillDiffusionMatrix(out);
        for (int i = 0; i < stateDimension; ++i) {
            for (int j = 0; j < stateDimension; ++j) {
                out[i][j] *= dt;
            }
        }
    }

    /**
     * Euler chain rule: ∂logL/∂A_{kl} += dt·(−(∂logL/∂F)_{kl} + (∂logL/∂f)_k·μ_l).
     */
    @Override
    public void accumulateSelectionGradient(final double dt,
                                            final double[][] dLogL_dF,
                                            final double[] dLogL_df,
                                            final double[] gradientAccumulator) {
        final double[] mu = new double[stateDimension];
        getInitialMean(mu);
        for (int k = 0; k < stateDimension; ++k) {
            for (int l = 0; l < stateDimension; ++l) {
                gradientAccumulator[k * stateDimension + l] +=
                        dt * (-dLogL_dF[k][l] + dLogL_df[k] * mu[l]);
            }
        }
    }

    @Override
    public void accumulateDiffusionGradient(final double dt,
                                            final double[][] dLogL_dV,
                                            final double[] gradientAccumulator) {
        for (int i = 0; i < stateDimension; ++i) {
            for (int j = 0; j < stateDimension; ++j) {
                gradientAccumulator[i * stateDimension + j] += dt * dLogL_dV[i][j];
            }
        }
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
    protected void storeState() { }

    @Override
    protected void restoreState() { }

    @Override
    protected void acceptState() { }

    private void validateShapes() {
        if (driftMatrix.getRowDimension() != stateDimension || driftMatrix.getColumnDimension() != stateDimension) {
            throw new IllegalArgumentException(
                    "driftMatrix must be stateDimension x stateDimension");
        }
        if (diffusionMatrix.getRowDimension() != stateDimension || diffusionMatrix.getColumnDimension() != stateDimension) {
            throw new IllegalArgumentException(
                    "diffusionMatrix must be stateDimension x stateDimension");
        }
        if (initialCovariance.getRowDimension() != stateDimension || initialCovariance.getColumnDimension() != stateDimension) {
            throw new IllegalArgumentException(
                    "initialCovariance must be stateDimension x stateDimension");
        }
        if (stationaryMean.getDimension() != 1 && stationaryMean.getDimension() != stateDimension) {
            throw new IllegalArgumentException(
                    "stationaryMean dimension must be 1 or stateDimension");
        }
    }

    private static double validatedDelta(final TimeGrid timeGrid,
                                         final int fromIndex,
                                         final int toIndex) {
        final double dt = timeGrid.getDelta(fromIndex, toIndex);
        if (!(dt > 0.0)) {
            throw new IllegalArgumentException("Time increments must be strictly positive");
        }
        return dt;
    }

    private static void checkVectorLength(final double[] vector,
                                          final int expectedLength,
                                          final String label) {
        if (vector == null || vector.length != expectedLength) {
            throw new IllegalArgumentException(label + " vector must have length " + expectedLength);
        }
    }

    private static void checkSquareMatrix(final double[][] matrix,
                                          final int expectedSize,
                                          final String label) {
        if (matrix == null || matrix.length != expectedSize) {
            throw new IllegalArgumentException(
                    label + " matrix must have size " + expectedSize + "x" + expectedSize);
        }
        for (double[] row : matrix) {
            if (row == null || row.length != expectedSize) {
                throw new IllegalArgumentException(
                        label + " matrix must have size " + expectedSize + "x" + expectedSize);
            }
        }
    }

    private static void setIdentity(final double[][] matrix) {
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = 0; j < matrix[i].length; ++j) {
                matrix[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }
    }

    private static void copyMatrixParameter(final MatrixParameterInterface parameter, final double[][] out) {
        for (int i = 0; i < out.length; ++i) {
            for (int j = 0; j < out[i].length; ++j) {
                out[i][j] = parameter.getParameterValue(i, j);
            }
        }
    }
}
