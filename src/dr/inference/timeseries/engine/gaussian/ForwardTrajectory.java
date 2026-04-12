package dr.inference.timeseries.engine.gaussian;

/**
 * Package-private snapshot of all quantities produced by the Kalman filter forward pass.
 *
 * <p>Transition arrays ({@code transitionMatrices}, {@code transitionOffsets},
 * {@code stepCovariances}) are indexed by the <em>from</em> step: index {@code t} holds
 * the parameters for the transition from time {@code t} to time {@code t+1}, so they
 * have length {@code timeCount − 1}.
 *
 * <p>Predicted and filtered arrays are indexed by absolute time step 0 … T−1.
 */
final class ForwardTrajectory {

    final int timeCount;
    final int stateDimension;

    /** m_{t|t}: filtered posterior mean at each step. */
    final double[][] filteredMeans;

    /** P_{t|t}: filtered posterior covariance at each step. */
    final double[][][] filteredCovariances;

    /**
     * m_{t|t−1}: one-step-ahead predicted mean at each step.
     * At t = 0 this equals the process initial mean.
     */
    final double[][] predictedMeans;

    /**
     * P_{t|t−1}: one-step-ahead predicted covariance at each step.
     * At t = 0 this equals the process initial covariance.
     */
    final double[][][] predictedCovariances;

    /**
     * F_t: transition matrix for the step from t to t+1.
     * Length: timeCount − 1.
     */
    final double[][][] transitionMatrices;

    /**
     * f_t: affine transition offset for the step from t to t+1.
     * Length: timeCount − 1.
     */
    final double[][] transitionOffsets;

    /**
     * Q_step_t = dt · Q: step noise covariance for the step from t to t+1.
     * Length: timeCount − 1.
     */
    final double[][][] stepCovariances;

    ForwardTrajectory(final int timeCount, final int stateDimension) {
        this.timeCount      = timeCount;
        this.stateDimension = stateDimension;
        final int steps = timeCount - 1;
        filteredMeans        = new double[timeCount][stateDimension];
        filteredCovariances  = new double[timeCount][stateDimension][stateDimension];
        predictedMeans       = new double[timeCount][stateDimension];
        predictedCovariances = new double[timeCount][stateDimension][stateDimension];
        transitionMatrices   = new double[steps][stateDimension][stateDimension];
        transitionOffsets    = new double[steps][stateDimension];
        stepCovariances      = new double[steps][stateDimension][stateDimension];
    }
}
