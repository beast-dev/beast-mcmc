package dr.inference.timeseries.engine.gaussian;

/**
 * Package-private container for the per-step sufficient statistics produced by the
 * Rauch-Tung-Striebel (RTS) backward pass.
 *
 * <p>The lag-1 smoothed cross-covariance is <em>not</em> stored directly. It is
 * recovered on demand from the identity
 * <pre>
 *   P_{t+1, t | T}  =  P_{t+1 | T} · G_t^T
 * </pre>
 * where {@code smootherGain} is G_t and {@code P_{t+1|T}} is the
 * {@code smoothedCovariance} of the <em>next</em> step's stats. This avoids an extra
 * d × d allocation per step while keeping gradient accumulation O(1) extra memory.
 */
final class BranchSmootherStats {

    /** Time index of this step. */
    final int timeIndex;

    /** m_{t|T} = E[x_t | Y_{1:T}]: smoothed posterior mean. */
    final double[] smoothedMean;

    /** P_{t|T} = Cov(x_t | Y_{1:T}): smoothed posterior covariance. */
    final double[][] smoothedCovariance;

    /**
     * G_t = P_{t|t} · F_t^T · P_{t+1|t}^{-1}: RTS smoother gain.
     *
     * <p>Null for the final time step (t = T − 1) since there is no forward transition.
     * Used to reconstruct P_{t+1,t|T} = P_{t+1|T} · G_t^T without storing it.
     */
    final double[][] smootherGain;

    BranchSmootherStats(final int timeIndex, final int stateDimension, final boolean hasForwardStep) {
        this.timeIndex      = timeIndex;
        smoothedMean        = new double[stateDimension];
        smoothedCovariance  = new double[stateDimension][stateDimension];
        smootherGain        = hasForwardStep ? new double[stateDimension][stateDimension] : null;
    }
}