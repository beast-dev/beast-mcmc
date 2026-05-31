package dr.inference.timeseries.engine.kalman;

import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;

/**
 * Snapshot of all quantities produced by the moment-form Kalman filter forward pass.
 *
 * <p>Transition arrays ({@code transitionMatrices}, {@code transitionOffsets},
 * {@code stepCovariances}) are indexed by the <em>from</em> step: index {@code t} holds
 * the parameters for the transition from time {@code t} to time {@code t+1}, so they
 * have length {@code timeCount − 1}.
 *
 * <p>Predicted and filtered arrays are indexed by absolute time step 0 … T−1.
 */
public final class ForwardTrajectory {

    public final int timeCount;
    public final int stateDimension;
    public final int matrixSize;

    /** m_{t|t}: filtered posterior mean at each step. */
    public final double[] filteredMeans;

    /** P_{t|t}: filtered posterior covariance at each step. */
    public final double[] filteredCovariances;

    /**
     * m_{t|t−1}: one-step-ahead predicted mean at each step.
     * At t = 0 this equals the process initial mean.
     */
    public final double[] predictedMeans;

    /**
     * P_{t|t−1}: one-step-ahead predicted covariance at each step.
     * At t = 0 this equals the process initial covariance.
     */
    public final double[] predictedCovariances;

    /**
     * F_t: transition matrix for the step from t to t+1.
     * Length: timeCount − 1.
     */
    public final double[] transitionMatrices;

    /**
     * f_t: affine transition offset for the step from t to t+1.
     * Length: timeCount − 1.
     */
    public final double[] transitionOffsets;

    /**
     * Q_step_t = dt · Q: step noise covariance for the step from t to t+1.
     * Length: timeCount − 1.
     */
    public final double[] stepCovariances;

    ForwardTrajectory(final int timeCount, final int stateDimension) {
        this.timeCount      = timeCount;
        this.stateDimension = stateDimension;
        this.matrixSize     = stateDimension * stateDimension;
        final int steps = timeCount - 1;
        filteredMeans        = new double[timeCount * stateDimension];
        filteredCovariances  = new double[timeCount * matrixSize];
        predictedMeans       = new double[timeCount * stateDimension];
        predictedCovariances = new double[timeCount * matrixSize];
        transitionMatrices   = new double[steps * matrixSize];
        transitionOffsets    = new double[steps * stateDimension];
        stepCovariances      = new double[steps * matrixSize];
    }

    public int stateVectorOffset(final int timeIndex) {
        checkTimeIndex(timeIndex);
        return timeIndex * stateDimension;
    }

    public int stateMatrixOffset(final int timeIndex) {
        checkTimeIndex(timeIndex);
        return timeIndex * matrixSize;
    }

    public int branchVectorOffset(final int branchIndex) {
        checkBranchIndex(branchIndex);
        return branchIndex * stateDimension;
    }

    public int branchMatrixOffset(final int branchIndex) {
        checkBranchIndex(branchIndex);
        return branchIndex * matrixSize;
    }

    void copyFilteredMeanFrom(final int timeIndex, final double[] source) {
        System.arraycopy(source, 0, filteredMeans, stateVectorOffset(timeIndex), stateDimension);
    }

    void copyFilteredCovarianceFrom(final int timeIndex, final double[] source) {
        System.arraycopy(source, 0, filteredCovariances, stateMatrixOffset(timeIndex), matrixSize);
    }

    void copyPredictedMeanFrom(final int timeIndex, final double[] source) {
        System.arraycopy(source, 0, predictedMeans, stateVectorOffset(timeIndex), stateDimension);
    }

    void copyPredictedCovarianceFrom(final int timeIndex, final double[] source) {
        System.arraycopy(source, 0, predictedCovariances, stateMatrixOffset(timeIndex), matrixSize);
    }

    void copyTransitionMatrixFrom(final int branchIndex, final double[] source) {
        System.arraycopy(source, 0, transitionMatrices, branchMatrixOffset(branchIndex), matrixSize);
    }

    void copyTransitionOffsetFrom(final int branchIndex, final double[] source) {
        System.arraycopy(source, 0, transitionOffsets, branchVectorOffset(branchIndex), stateDimension);
    }

    void copyStepCovarianceFrom(final int branchIndex, final double[] source) {
        System.arraycopy(source, 0, stepCovariances, branchMatrixOffset(branchIndex), matrixSize);
    }

    public void copyFilteredMeanTo(final int timeIndex, final double[] out) {
        System.arraycopy(filteredMeans, stateVectorOffset(timeIndex), out, 0, stateDimension);
    }

    public void copyFilteredCovarianceTo(final int timeIndex, final double[][] out) {
        MatrixOps.fromFlat(filteredCovariances, stateMatrixOffset(timeIndex), out, stateDimension);
    }

    public void copyPredictedMeanTo(final int timeIndex, final double[] out) {
        System.arraycopy(predictedMeans, stateVectorOffset(timeIndex), out, 0, stateDimension);
    }

    public void copyPredictedCovarianceTo(final int timeIndex, final double[][] out) {
        MatrixOps.fromFlat(predictedCovariances, stateMatrixOffset(timeIndex), out, stateDimension);
    }

    public void copyTransitionMatrixTo(final int branchIndex, final double[][] out) {
        MatrixOps.fromFlat(transitionMatrices, branchMatrixOffset(branchIndex), out, stateDimension);
    }

    public void copyTransitionOffsetTo(final int branchIndex, final double[] out) {
        System.arraycopy(transitionOffsets, branchVectorOffset(branchIndex), out, 0, stateDimension);
    }

    public void copyStepCovarianceTo(final int branchIndex, final double[][] out) {
        MatrixOps.fromFlat(stepCovariances, branchMatrixOffset(branchIndex), out, stateDimension);
    }

    private void checkTimeIndex(final int timeIndex) {
        if (timeIndex < 0 || timeIndex >= timeCount) {
            throw new IllegalArgumentException("timeIndex out of bounds: " + timeIndex);
        }
    }

    private void checkBranchIndex(final int branchIndex) {
        if (branchIndex < 0 || branchIndex >= timeCount - 1) {
            throw new IllegalArgumentException("branchIndex out of bounds: " + branchIndex);
        }
    }
}
