package dr.inference.timeseries.representation;

/**
 * Internal algebra used to represent Gaussian quantities.
 *
 * <p>{@link #EXPECTATION} corresponds to the familiar mean/covariance form used by
 * the current Kalman code. {@link #CANONICAL} corresponds to precision/information
 * form and is intended for future high-performance pathways that avoid repeated
 * expectation-form conversions inside inference engines.
 */
public enum GaussianComputationMode {
    EXPECTATION,
    CANONICAL
}
