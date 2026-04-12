package dr.inference.timeseries.engine.gaussian;

import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.engine.LikelihoodEngine;
import dr.inference.timeseries.gaussian.GaussianObservationModel;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;

/**
 * Linear-Gaussian likelihood engine based on a Kalman filter.
 *
 * <p>This engine marginalizes over the latent state trajectory. Missing observations
 * are supported by encoding an entire observation column as {@link Double#NaN}.
 *
 * <p>The forward pass is factored into {@link #runForwardPass(ForwardTrajectory)} so
 * that subclasses (e.g. {@link KalmanSmootherEngine}) can store the full trajectory
 * for a subsequent RTS backward pass without duplicating the filter logic. When
 * {@code trajectoryOut} is {@code null} the trajectory is computed but not stored,
 * keeping the base-class path allocation-free.
 *
 * <p>All static linear-algebra helpers are package-private so that other classes in
 * this package (smoother, gradient formulas) can reuse them without duplication.
 */
public class KalmanLikelihoodEngine implements LikelihoodEngine {

    static final double LOG_TWO_PI = Math.log(2.0 * Math.PI);
    private static final double SYMMETRY_TOLERANCE = 1E-12;
    static final double MIN_DIAGONAL_JITTER = 1E-10;

    private final GaussianTransitionRepresentation transitionRepresentation;
    private final GaussianObservationModel observationModel;
    private final TimeGrid timeGrid;

    final int stateDimension;
    final int observationDimension;
    final int timeCount;

    private final double[] filteredMean;
    private final double[] predictedMean;
    private final double[] offset;
    private final double[] innovation;
    private final double[] kalmanGainColumn;
    private final double[] workingVector;
    private final double[] observationVector;

    private final double[][] filteredCovariance;
    private final double[][] predictedCovariance;
    private final double[][] transitionMatrix;
    private final double[][] transitionCovariance;
    private final double[][] initialCovariance;
    private final double[][] designMatrix;
    private final double[][] noiseCovariance;
    private final double[][] innovationCovariance;
    private final double[][] innovationCovarianceInverse;
    private final double[][] gainMatrix;
    private final double[][] tempStateState;
    private final double[][] tempStateObs;
    private final double[][] tempObsState;
    private final double[][] tempObsObs;

    private boolean likelihoodKnown;
    private double logLikelihood;

    public KalmanLikelihoodEngine(final GaussianTransitionRepresentation transitionRepresentation,
                                  final GaussianObservationModel observationModel,
                                  final TimeGrid timeGrid) {
        if (transitionRepresentation == null) {
            throw new IllegalArgumentException("transitionRepresentation must not be null");
        }
        if (observationModel == null) {
            throw new IllegalArgumentException("observationModel must not be null");
        }
        if (timeGrid == null) {
            throw new IllegalArgumentException("timeGrid must not be null");
        }
        this.transitionRepresentation = transitionRepresentation;
        this.observationModel         = observationModel;
        this.timeGrid                 = timeGrid;
        this.stateDimension           = transitionRepresentation.getStateDimension();
        this.observationDimension     = observationModel.getObservationDimension();
        this.timeCount                = timeGrid.getTimeCount();

        if (observationModel.getTimeCount() != timeCount) {
            throw new IllegalArgumentException(
                    "Observation time count (" + observationModel.getTimeCount() +
                    ") must match time grid count (" + timeCount + ")");
        }

        filteredMean              = new double[stateDimension];
        predictedMean             = new double[stateDimension];
        offset                    = new double[stateDimension];
        innovation                = new double[observationDimension];
        kalmanGainColumn          = new double[stateDimension];
        workingVector             = new double[stateDimension];
        observationVector         = new double[observationDimension];

        filteredCovariance        = new double[stateDimension][stateDimension];
        predictedCovariance       = new double[stateDimension][stateDimension];
        transitionMatrix          = new double[stateDimension][stateDimension];
        transitionCovariance      = new double[stateDimension][stateDimension];
        initialCovariance         = new double[stateDimension][stateDimension];
        designMatrix              = new double[observationDimension][stateDimension];
        noiseCovariance           = new double[observationDimension][observationDimension];
        innovationCovariance      = new double[observationDimension][observationDimension];
        innovationCovarianceInverse = new double[observationDimension][observationDimension];
        gainMatrix                = new double[stateDimension][observationDimension];
        tempStateState            = new double[stateDimension][stateDimension];
        tempStateObs              = new double[stateDimension][observationDimension];
        tempObsState              = new double[observationDimension][stateDimension];
        tempObsObs                = new double[observationDimension][observationDimension];
    }

    // ─── LikelihoodEngine ───────────────────────────────────────────────────────

    @Override
    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = computeLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
    }

    // ─── Accessors ──────────────────────────────────────────────────────────────

    public GaussianTransitionRepresentation getTransitionRepresentation() {
        return transitionRepresentation;
    }

    public GaussianObservationModel getObservationModel() {
        return observationModel;
    }

    public TimeGrid getTimeGrid() {
        return timeGrid;
    }

    // ─── Forward pass ───────────────────────────────────────────────────────────

    /**
     * Hook for subclasses: by default delegates to {@link #runForwardPass(ForwardTrajectory)}
     * with a null trajectory (no storage overhead).
     */
    protected double computeLogLikelihood() {
        return runForwardPass(null);
    }

    /**
     * Runs the Kalman filter forward pass and returns the log-likelihood.
     *
     * <p>When {@code trajectoryOut} is non-null, filtered/predicted moments and transition
     * parameters are written into it at every step. This is the only place where trajectory
     * storage is performed; the rest of the engine is unaffected.
     *
     * @param trajectoryOut destination for the full trajectory, or {@code null} to skip storage
     * @return log p(Y | θ) accumulated over all observed steps
     */
    protected double runForwardPass(final ForwardTrajectory trajectoryOut) {
        observationModel.fillDesignMatrix(designMatrix);
        observationModel.fillNoiseCovariance(noiseCovariance);
        transitionRepresentation.getInitialMean(filteredMean);
        transitionRepresentation.getInitialCovariance(initialCovariance);
        copyMatrix(initialCovariance, filteredCovariance);

        double value = 0.0;

        for (int timeIndex = 0; timeIndex < timeCount; ++timeIndex) {

            // ── Prediction ──────────────────────────────────────────────────────
            if (timeIndex == 0) {
                copyVector(filteredMean, predictedMean);
                copyMatrix(filteredCovariance, predictedCovariance);
            } else {
                transitionRepresentation.getTransitionMatrix(
                        timeIndex - 1, timeIndex, timeGrid, transitionMatrix);
                transitionRepresentation.getTransitionOffset(
                        timeIndex - 1, timeIndex, timeGrid, offset);
                transitionRepresentation.getTransitionCovariance(
                        timeIndex - 1, timeIndex, timeGrid, transitionCovariance);

                multiplyMatrixVector(transitionMatrix, filteredMean, predictedMean);
                addInPlace(predictedMean, offset);

                multiplyMatrixMatrix(transitionMatrix, filteredCovariance, tempStateState);
                multiplyMatrixMatrixTransposedRight(tempStateState, transitionMatrix, predictedCovariance);
                addMatrixInPlace(predictedCovariance, transitionCovariance);
                symmetrize(predictedCovariance);

                if (trajectoryOut != null) {
                    copyMatrix(transitionMatrix,    trajectoryOut.transitionMatrices[timeIndex - 1]);
                    copyVector(offset,              trajectoryOut.transitionOffsets[timeIndex - 1]);
                    copyMatrix(transitionCovariance, trajectoryOut.stepCovariances[timeIndex - 1]);
                }
            }

            if (trajectoryOut != null) {
                copyVector(predictedMean,       trajectoryOut.predictedMeans[timeIndex]);
                copyMatrix(predictedCovariance, trajectoryOut.predictedCovariances[timeIndex]);
            }

            // ── Update ──────────────────────────────────────────────────────────
            if (observationModel.isObservationMissing(timeIndex)) {
                copyVector(predictedMean,       filteredMean);
                copyMatrix(predictedCovariance, filteredCovariance);
            } else {
                observationModel.fillObservationVector(timeIndex, observationVector);

                multiplyMatrixVector(designMatrix, predictedMean, innovation);
                for (int i = 0; i < observationDimension; ++i) {
                    innovation[i] = observationVector[i] - innovation[i];
                }

                multiplyMatrixMatrix(designMatrix, predictedCovariance, tempObsState);
                multiplyMatrixMatrixTransposedRight(tempObsState, designMatrix, innovationCovariance);
                addMatrixInPlace(innovationCovariance, noiseCovariance);
                symmetrize(innovationCovariance);
                addJitterToDiagonal(innovationCovariance, MIN_DIAGONAL_JITTER);

                final CholeskyFactor innovationCholesky = cholesky(innovationCovariance);
                final double logDetInnovation = innovationCholesky.logDeterminant();
                copyMatrix(innovationCovariance, innovationCovarianceInverse);
                invertPositiveDefiniteFromCholesky(innovationCovarianceInverse, innovationCholesky);

                value += -0.5 * (observationDimension * LOG_TWO_PI
                        + logDetInnovation
                        + quadraticForm(innovationCovarianceInverse, innovation));

                multiplyMatrixMatrixTransposedRight(predictedCovariance, designMatrix, tempStateObs);
                multiplyMatrixMatrix(tempStateObs, innovationCovarianceInverse, gainMatrix,
                        stateDimension, observationDimension, observationDimension);

                multiplyMatrixVector(gainMatrix, innovation, workingVector, stateDimension, observationDimension);
                for (int i = 0; i < stateDimension; ++i) {
                    filteredMean[i] = predictedMean[i] + workingVector[i];
                }

                // Joseph-stabilised covariance update: P = (I−KH) P (I−KH)^T + K R K^T
                //
                // predictedCovariance is no longer needed after the gain computation above,
                // so it is reused as a [d×d] scratch buffer. This avoids the shape mismatch
                // that would arise from using tempObsObs ([p×p]) for a [d×d] intermediate.
                multiplyMatrixMatrix(gainMatrix, designMatrix, tempStateState,
                        stateDimension, observationDimension, stateDimension);
                identityMinus(tempStateState);
                multiplyMatrixMatrix(tempStateState, predictedCovariance, filteredCovariance,
                        stateDimension, stateDimension, stateDimension);
                multiplyMatrixMatrixTransposedRight(filteredCovariance, tempStateState, predictedCovariance);
                copyMatrix(predictedCovariance, filteredCovariance);
                multiplyMatrixMatrix(gainMatrix, noiseCovariance, tempStateObs,
                        stateDimension, observationDimension, observationDimension);
                multiplyMatrixMatrixTransposedRight(tempStateObs, gainMatrix, tempStateState);
                addMatrixInPlace(filteredCovariance, tempStateState);
                symmetrize(filteredCovariance);
            }

            if (trajectoryOut != null) {
                copyVector(filteredMean,       trajectoryOut.filteredMeans[timeIndex]);
                copyMatrix(filteredCovariance, trajectoryOut.filteredCovariances[timeIndex]);
            }
        }

        return value;
    }

    // ─── Package-private linear-algebra utilities ────────────────────────────────
    // Package-private so that KalmanSmootherEngine and gradient formulas in this
    // package can share them without duplication.

    static void multiplyMatrixVector(final double[][] matrix,
                                     final double[] vector,
                                     final double[] out) {
        multiplyMatrixVector(matrix, vector, out, matrix.length, vector.length);
    }

    static void multiplyMatrixVector(final double[][] matrix,
                                     final double[] vector,
                                     final double[] out,
                                     final int rows,
                                     final int cols) {
        for (int i = 0; i < rows; ++i) {
            double sum = 0.0;
            for (int j = 0; j < cols; ++j) {
                sum += matrix[i][j] * vector[j];
            }
            out[i] = sum;
        }
    }

    static void multiplyMatrixMatrix(final double[][] left,
                                     final double[][] right,
                                     final double[][] out) {
        multiplyMatrixMatrix(left, right, out, left.length, right.length, right[0].length);
    }

    static void multiplyMatrixMatrix(final double[][] left,
                                     final double[][] right,
                                     final double[][] out,
                                     final int rows,
                                     final int inner,
                                     final int cols) {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double sum = 0.0;
                for (int k = 0; k < inner; ++k) {
                    sum += left[i][k] * right[k][j];
                }
                out[i][j] = sum;
            }
        }
    }

    /** Computes {@code out = left · right^T}. */
    static void multiplyMatrixMatrixTransposedRight(final double[][] left,
                                                    final double[][] right,
                                                    final double[][] out) {
        final int rows  = left.length;
        final int inner = left[0].length;
        final int cols  = right.length;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double sum = 0.0;
                for (int k = 0; k < inner; ++k) {
                    sum += left[i][k] * right[j][k];
                }
                out[i][j] = sum;
            }
        }
    }

    static void addMatrixInPlace(final double[][] accumulator, final double[][] increment) {
        for (int i = 0; i < accumulator.length; ++i) {
            for (int j = 0; j < accumulator[i].length; ++j) {
                accumulator[i][j] += increment[i][j];
            }
        }
    }

    static void addInPlace(final double[] accumulator, final double[] increment) {
        for (int i = 0; i < accumulator.length; ++i) {
            accumulator[i] += increment[i];
        }
    }

    static void identityMinus(final double[][] matrix) {
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = 0; j < matrix[i].length; ++j) {
                matrix[i][j] = (i == j ? 1.0 : 0.0) - matrix[i][j];
            }
        }
    }

    static void copyVector(final double[] source, final double[] target) {
        System.arraycopy(source, 0, target, 0, source.length);
    }

    static void copyMatrix(final double[][] source, final double[][] target) {
        copyMatrix(source, target, source.length, source[0].length);
    }

    static void copyMatrix(final double[][] source,
                           final double[][] target,
                           final int rows,
                           final int cols) {
        for (int i = 0; i < rows; ++i) {
            System.arraycopy(source[i], 0, target[i], 0, cols);
        }
    }

    static void symmetrize(final double[][] matrix) {
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = i + 1; j < matrix[i].length; ++j) {
                final double average = 0.5 * (matrix[i][j] + matrix[j][i]);
                matrix[i][j] = average;
                matrix[j][i] = average;
            }
        }
    }

    static void addJitterToDiagonal(final double[][] matrix, final double minimumJitter) {
        for (int i = 0; i < matrix.length; ++i) {
            if (matrix[i][i] < minimumJitter) {
                matrix[i][i] = minimumJitter;
            }
        }
    }

    static double quadraticForm(final double[][] matrix, final double[] vector) {
        double result = 0.0;
        for (int i = 0; i < vector.length; ++i) {
            double rowSum = 0.0;
            for (int j = 0; j < vector.length; ++j) {
                rowSum += matrix[i][j] * vector[j];
            }
            result += vector[i] * rowSum;
        }
        return result;
    }

    static CholeskyFactor cholesky(final double[][] matrix) {
        final int dim = matrix.length;
        final double[][] lower = new double[dim][dim];

        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j <= i; ++j) {
                double sum = matrix[i][j];
                for (int k = 0; k < j; ++k) {
                    sum -= lower[i][k] * lower[j][k];
                }

                if (i == j) {
                    if (sum <= 0.0) {
                        throw new IllegalArgumentException("Matrix is not positive definite");
                    }
                    lower[i][j] = Math.sqrt(sum);
                } else {
                    lower[i][j] = sum / lower[j][j];
                }
            }

            for (int j = i + 1; j < dim; ++j) {
                if (Math.abs(matrix[i][j] - matrix[j][i]) > SYMMETRY_TOLERANCE) {
                    throw new IllegalArgumentException("Matrix must be symmetric");
                }
            }
        }

        return new CholeskyFactor(lower);
    }

    static void invertPositiveDefiniteFromCholesky(final double[][] out,
                                                   final CholeskyFactor factor) {
        final int dim = out.length;
        final double[][] inverse = new double[dim][dim];
        final double[] basis    = new double[dim];
        final double[] solution = new double[dim];

        for (int column = 0; column < dim; ++column) {
            for (int i = 0; i < dim; ++i) {
                basis[i] = 0.0;
            }
            basis[column] = 1.0;
            factor.solveSymmetricSystem(basis, solution);
            for (int row = 0; row < dim; ++row) {
                inverse[row][column] = solution[row];
            }
        }

        copyMatrix(inverse, out);
        symmetrize(out);
    }

    // ─── Cholesky factor ────────────────────────────────────────────────────────

    static final class CholeskyFactor {
        private final double[][] lower;

        CholeskyFactor(final double[][] lower) {
            this.lower = lower;
        }

        double logDeterminant() {
            double value = 0.0;
            for (int i = 0; i < lower.length; ++i) {
                value += 2.0 * Math.log(lower[i][i]);
            }
            return value;
        }

        void solveSymmetricSystem(final double[] rhs, final double[] out) {
            final int dim = lower.length;
            final double[] y = new double[dim];

            for (int i = 0; i < dim; ++i) {
                double sum = rhs[i];
                for (int j = 0; j < i; ++j) {
                    sum -= lower[i][j] * y[j];
                }
                y[i] = sum / lower[i][i];
            }

            for (int i = dim - 1; i >= 0; --i) {
                double sum = y[i];
                for (int j = i + 1; j < dim; ++j) {
                    sum -= lower[j][i] * out[j];
                }
                out[i] = sum / lower[i][i];
            }
        }
    }
}
