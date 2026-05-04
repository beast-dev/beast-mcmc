package dr.inference.timeseries.engine.gaussian;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;

/**
 * Canonical-form snapshot of the forward and backward Gaussian trajectory.
 *
 * <p>This structure intentionally omits the expectation-form filtered/predicted
 * means and covariances used by the legacy Kalman pathway. It is the native
 * trajectory surface for canonical smoothing and canonical analytical gradients.
 */
public final class CanonicalForwardTrajectory {

    public final int timeCount;
    public final int stateDimension;

    public final CanonicalGaussianState[] filteredStates;
    public final CanonicalGaussianState[] predictedStates;
    public final CanonicalGaussianState[] smoothedStates;
    public final CanonicalGaussianTransition[] transitions;
    public final CanonicalGaussianState[] branchPairStates;

    CanonicalForwardTrajectory(final int timeCount, final int stateDimension) {
        this.timeCount = timeCount;
        this.stateDimension = stateDimension;
        final int steps = timeCount - 1;
        filteredStates = new CanonicalGaussianState[timeCount];
        predictedStates = new CanonicalGaussianState[timeCount];
        smoothedStates = new CanonicalGaussianState[timeCount];
        transitions = new CanonicalGaussianTransition[steps];
        branchPairStates = new CanonicalGaussianState[steps];
        for (int t = 0; t < timeCount; ++t) {
            filteredStates[t] = new CanonicalGaussianState(stateDimension);
            predictedStates[t] = new CanonicalGaussianState(stateDimension);
            smoothedStates[t] = new CanonicalGaussianState(stateDimension);
            if (t < steps) {
                transitions[t] = new CanonicalGaussianTransition(stateDimension);
                branchPairStates[t] = new CanonicalGaussianState(2 * stateDimension);
            }
        }
    }
}
