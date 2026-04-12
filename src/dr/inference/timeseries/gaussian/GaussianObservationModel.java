package dr.inference.timeseries.gaussian;

import dr.inference.model.AbstractModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.timeseries.core.ObservationModel;

/**
 * Gaussian observation model for a linear state-space system.
 *
 * Observations are stored in a matrix whose columns correspond to time points and whose rows
 * correspond to observed dimensions. Missing values may be encoded with {@link Double#NaN}.
 */
public class GaussianObservationModel extends AbstractModel implements ObservationModel {

    private final int observationDimension;
    private final MatrixParameter designMatrix;
    private final MatrixParameter noiseCovariance;
    private final MatrixParameter observations;
    private final int timeCount;

    public GaussianObservationModel(final String name,
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
        this.observationDimension = observationDimension;
        this.designMatrix = designMatrix;
        this.noiseCovariance = noiseCovariance;
        this.observations = observations;
        this.timeCount = inferTimeCount(observations, observationDimension);

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

    public void fillDesignMatrix(final double[][] out) {
        checkRectangularMatrix(out, observationDimension, out[0].length, "designMatrix output");
        for (int i = 0; i < out.length; ++i) {
            for (int j = 0; j < out[i].length; ++j) {
                out[i][j] = designMatrix.getParameterValue(i, j);
            }
        }
    }

    public void fillNoiseCovariance(final double[][] out) {
        checkRectangularMatrix(out, observationDimension, observationDimension, "noise covariance output");
        for (int i = 0; i < observationDimension; ++i) {
            for (int j = 0; j < observationDimension; ++j) {
                out[i][j] = noiseCovariance.getParameterValue(i, j);
            }
        }
    }

    public void fillObservationVector(final int timeIndex, final double[] out) {
        if (timeIndex < 0 || timeIndex >= timeCount) {
            throw new IllegalArgumentException("timeIndex out of bounds: " + timeIndex);
        }
        if (out == null || out.length != observationDimension) {
            throw new IllegalArgumentException("observation vector must have length " + observationDimension);
        }

        final int nonMissingCount = countNonMissingEntries(timeIndex);
        if (nonMissingCount > 0 && nonMissingCount < observationDimension) {
            throw new IllegalArgumentException(
                    "Partial missing observations are not currently supported at time index " + timeIndex +
                    "; either provide a complete observation vector or mark the whole column as NaN.");
        }

        for (int i = 0; i < observationDimension; ++i) {
            out[i] = observations.getParameterValue(i, timeIndex);
        }
    }

    public boolean isObservationMissing(final int timeIndex) {
        if (timeIndex < 0 || timeIndex >= timeCount) {
            throw new IllegalArgumentException("timeIndex out of bounds: " + timeIndex);
        }
        return countNonMissingEntries(timeIndex) == 0;
    }


    private int countNonMissingEntries(final int timeIndex) {
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

    private static int inferTimeCount(final MatrixParameter observations, final int observationDimension) {
        final int totalDimension = observations.getDimension();
        if (totalDimension % observationDimension != 0) {
            throw new IllegalArgumentException(
                    "Observation matrix flattened dimension (" + totalDimension +
                    ") must be divisible by observationDimension (" + observationDimension + ")");
        }
        return totalDimension / observationDimension;
    }

    private static void checkRectangularMatrix(final double[][] matrix,
                                               final int expectedRows,
                                               final int expectedCols,
                                               final String label) {
        if (matrix == null || matrix.length != expectedRows) {
            throw new IllegalArgumentException(label + " must have " + expectedRows + " rows");
        }
        for (double[] row : matrix) {
            if (row == null || row.length != expectedCols) {
                throw new IllegalArgumentException(label + " must have " + expectedCols + " columns");
            }
        }
    }
}
