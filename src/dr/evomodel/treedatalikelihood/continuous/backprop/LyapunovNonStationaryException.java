package dr.evomodel.treedatalikelihood.continuous.backprop;

/**
 * Signals that the continuous Lyapunov system is singular or non-stationary for the
 * current block-diagonal drift parametrization.
 */
public final class LyapunovNonStationaryException extends RuntimeException {

    public LyapunovNonStationaryException(final String message) {
        super(message);
    }
}
