package dr.inference.timeseries.model.gaussian;

import dr.inference.model.AbstractModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.timeseries.core.ObservationModel;

/**
 * Gaussian observation model for a linear state-space system,
 * {@code y_t = H x_t + e_t} with {@code e_t ~ N(0, R)}.
 *
 * <p><b>Layout.</b> All three matrices follow the {@link MatrixParameter} convention that each
 * component parameter is one column:
 * <ul>
 *   <li>{@code designMatrix} ({@code H}) is {@code observationDimension x stateDimension};</li>
 *   <li>{@code noiseCovariance} ({@code R}) is {@code observationDimension x observationDimension};</li>
 *   <li>{@code observations} ({@code Y}) is {@code observationDimension x timeCount}: column
 *       {@code t} holds the observation vector at time index {@code t}.  In XML each
 *       {@code <parameter>} child of the observation {@code <matrixParameter>} is therefore one time
 *       point of length {@code observationDimension}, not one observed coordinate across time.</li>
 * </ul>
 * The constructor rejects every other shape with a message naming the offending matrix.  The state
 * dimension is not known to this class; callers that know the latent process check the design
 * matrix column count through {@link #validateStateDimension(int)}.
 *
 * <p><b>Missingness.</b> A missing value is encoded as {@link Double#NaN} in {@code observations}.
 * The canonical engines, which are the production forward, smoother and gradient modes, read the
 * observation matrix directly and handle an arbitrary per-cell NaN mask through
 * {@code CanonicalObservationUpdate}: any subset of the coordinates of a column may be observed,
 * including none.  The moment-form debug engines ({@code MomentKalmanLikelihoodEngine} and the
 * smoother and gradient engines built on its forward pass) instead go through
 * {@link #isObservationColumnMissing(int)} and {@link #fillObservationVector(int, double[])},
 * which distinguish only fully observed from fully missing columns; a partially observed column
 * is rejected there with an {@link IllegalArgumentException} that points to the canonical modes.
 */
public class LinearGaussianObservationModel extends AbstractModel implements ObservationModel {

    private final int observationDimension;
    private final MatrixParameter designMatrix;
    private final MatrixParameter noiseCovariance;
    private final MatrixParameter observations;
    private final int timeCount;

    public LinearGaussianObservationModel(final String name,
                                          final int observationDimension,
                                          final MatrixParameter designMatrix,
                                          final MatrixParameter noiseCovariance,
                                          final MatrixParameter observations) {
        super(name);
        if (observationDimension < 1) {
            throw new IllegalArgumentException("observationDimension must be at least 1");
        }
        if (designMatrix == null) {
            throw new IllegalArgumentException("designMatrix must not be null");
        }
        if (noiseCovariance == null) {
            throw new IllegalArgumentException("noiseCovariance must not be null");
        }
        if (observations == null) {
            throw new IllegalArgumentException("observations must not be null");
        }
        validateDesignMatrixShape(designMatrix, observationDimension);
        validateNoiseCovarianceShape(noiseCovariance, observationDimension);
        this.timeCount = validateObservationsShape(observations, observationDimension);

        this.observationDimension = observationDimension;
        this.designMatrix = designMatrix;
        this.noiseCovariance = noiseCovariance;
        this.observations = observations;

        addVariable(designMatrix);
        addVariable(noiseCovariance);
        addVariable(observations);
    }

    @Override
    public int getObservationDimension() {
        return observationDimension;
    }

    public int getTimeCount() {
        return timeCount;
    }

    public MatrixParameter getDesignMatrix() {
        return designMatrix;
    }

    public MatrixParameter getNoiseCovariance() {
        return noiseCovariance;
    }

    public MatrixParameter getObservations() {
        return observations;
    }

    /**
     * Checks that the design matrix has exactly one column per latent state coordinate.
     *
     * <p>The observation model alone does not know the latent state dimension, so this check is
     * run where the latent process is available: the {@code timeSeriesModel} parser and the
     * Gaussian time-series likelihood factory.
     *
     * @throws IllegalArgumentException when the design matrix column count differs from
     *                                  {@code stateDimension}
     */
    public void validateStateDimension(final int stateDimension) {
        if (stateDimension < 1) {
            throw new IllegalArgumentException("stateDimension must be at least 1");
        }
        final int columns = designMatrix.getColumnDimension();
        if (columns != stateDimension) {
            throw new IllegalArgumentException(
                    "designMatrix must be observationDimension x stateDimension = "
                            + observationDimension + " x " + stateDimension
                            + " (one column per latent state coordinate); found "
                            + describeShape(designMatrix) + ".");
        }
    }

    public void fillDesignMatrixFlat(final double[] out,
                                     final int stateDimension) {
        checkFlatMatrix(out, observationDimension, stateDimension, "designMatrix output");
        for (int i = 0; i < observationDimension; ++i) {
            final int rowOffset = i * stateDimension;
            for (int j = 0; j < stateDimension; ++j) {
                out[rowOffset + j] = designMatrix.getParameterValue(i, j);
            }
        }
    }

    public void fillNoiseCovarianceFlat(final double[] out) {
        checkFlatMatrix(out, observationDimension, observationDimension, "noise covariance output");
        for (int i = 0; i < observationDimension; ++i) {
            final int rowOffset = i * observationDimension;
            for (int j = 0; j < observationDimension; ++j) {
                out[rowOffset + j] = noiseCovariance.getParameterValue(i, j);
            }
        }
    }

    /**
     * Copies the observation vector at {@code timeIndex} into {@code out}.
     *
     * <p><b>Complete-column contract.</b> This accessor serves the moment-form debug engines,
     * which cannot mask individual coordinates.  The column must be fully observed or fully
     * missing; the usual call pattern tests {@link #isObservationColumnMissing(int)} first and
     * skips the measurement update when it returns {@code true}.  A partially observed column is
     * rejected.  The canonical engines never call this method: they read
     * {@link #getObservations()} directly and apply a per-cell NaN mask, so partially observed
     * columns are fully supported in the canonical forward, smoother and gradient modes.
     *
     * @throws IllegalArgumentException when {@code timeIndex} is out of range, {@code out} has the
     *                                  wrong length, or the column is partially observed
     */
    public void fillObservationVector(final int timeIndex, final double[] out) {
        checkTimeIndex(timeIndex);
        if (out == null || out.length != observationDimension) {
            throw new IllegalArgumentException("observation vector must have length " + observationDimension);
        }

        final int observedCount = countObservedEntries(timeIndex);
        if (observedCount > 0 && observedCount < observationDimension) {
            throw new IllegalArgumentException(
                    "Observation column " + timeIndex + " is partially observed (" + observedCount
                            + " of " + observationDimension + " coordinates). Partially observed columns are"
                            + " supported only by the canonical forward, smoother and gradient modes"
                            + " (CanonicalObservationUpdate masks NaN cells per column); the moment-form"
                            + " debug engines require every column to be fully observed or fully NaN."
                            + " Use the canonical modes (the timeSeriesLikelihood defaults) or mark the"
                            + " whole column NaN.");
        }

        for (int i = 0; i < observationDimension; ++i) {
            out[i] = observations.getParameterValue(i, timeIndex);
        }
    }

    /**
     * @return {@code true} when every coordinate of observation column {@code timeIndex} is NaN.
     *         A column with some but not all coordinates observed returns {@code false}; such a
     *         column is handled by the canonical engines and rejected by
     *         {@link #fillObservationVector(int, double[])}.
     */
    public boolean isObservationColumnMissing(final int timeIndex) {
        checkTimeIndex(timeIndex);
        return countObservedEntries(timeIndex) == 0;
    }

    private void checkTimeIndex(final int timeIndex) {
        if (timeIndex < 0 || timeIndex >= timeCount) {
            throw new IllegalArgumentException("timeIndex out of bounds: " + timeIndex);
        }
    }

    private int countObservedEntries(final int timeIndex) {
        int count = 0;
        for (int i = 0; i < observationDimension; ++i) {
            if (!Double.isNaN(observations.getParameterValue(i, timeIndex))) {
                count++;
            }
        }
        return count;
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

    // ------------------------------------------------------------------ shape validation

    private static void validateDesignMatrixShape(final MatrixParameter designMatrix,
                                                  final int observationDimension) {
        if (designMatrix.getColumnDimension() < 1
                || designMatrix.getRowDimension() != observationDimension) {
            throw new IllegalArgumentException(
                    "designMatrix must have observationDimension = " + observationDimension
                            + " rows (one per observed coordinate) and one column per latent state"
                            + " coordinate; found " + describeShape(designMatrix) + ". In XML each"
                            + " <parameter> child of the design <matrixParameter> is one column of H"
                            + " with length observationDimension.");
        }
    }

    private static void validateNoiseCovarianceShape(final MatrixParameter noiseCovariance,
                                                     final int observationDimension) {
        if (noiseCovariance.getColumnDimension() != observationDimension
                || noiseCovariance.getRowDimension() != observationDimension) {
            throw new IllegalArgumentException(
                    "noiseCovariance must be square with observationDimension = " + observationDimension
                            + " rows and columns; found " + describeShape(noiseCovariance) + ".");
        }
    }

    /**
     * @return the number of time points (columns) of a correctly shaped observation matrix
     */
    private static int validateObservationsShape(final MatrixParameter observations,
                                                 final int observationDimension) {
        final int columns = observations.getColumnDimension();
        if (columns < 1) {
            throw new IllegalArgumentException("observations must contain at least one time point (column)");
        }
        final int rows = observations.getRowDimension();
        if (rows != observationDimension) {
            final StringBuilder message = new StringBuilder()
                    .append("observations must have observationDimension = ").append(observationDimension)
                    .append(" rows and one column per time point; found ")
                    .append(describeShape(observations))
                    .append(". In XML each <parameter> child of the observation <matrixParameter> is one")
                    .append(" time point (a column) holding the ").append(observationDimension)
                    .append(" observed coordinates at that time.");
            if (columns == observationDimension) {
                message.append(" The supplied matrix has observationDimension columns and ").append(rows)
                        .append(" rows, so it looks transposed: each <parameter> child currently spans")
                        .append(" the time axis instead of the observed coordinates.");
            }
            throw new IllegalArgumentException(message.toString());
        }
        if (observations.getDimension() != rows * columns) {
            throw new IllegalArgumentException(
                    "observations must be rectangular: expected " + rows + " x " + columns + " = "
                            + (rows * columns) + " entries but the matrix holds "
                            + observations.getDimension() + ".");
        }
        return columns;
    }

    private static String describeShape(final MatrixParameter matrix) {
        final int columns = matrix.getColumnDimension();
        final String rows = columns < 1 ? "0" : String.valueOf(matrix.getRowDimension());
        return rows + " x " + columns + " (rows x columns)";
    }

    private static void checkFlatMatrix(final double[] matrix,
                                        final int expectedRows,
                                        final int expectedCols,
                                        final String label) {
        if (matrix == null || matrix.length != expectedRows * expectedCols) {
            throw new IllegalArgumentException(
                    label + " must have length " + (expectedRows * expectedCols));
        }
    }
}
