package dr.inference.timeseries.representation;

/**
 * Internal algebra used to represent Gaussian quantities.
 *
 * <p>{@link #MOMENT} corresponds to the familiar mean/covariance form used by
 * the current Kalman code. {@link #CANONICAL} corresponds to precision/information
 * form and is intended for future high-performance pathways that avoid repeated
 * moment-form conversions inside inference engines.
 */
public enum GaussianComputationMode {
    MOMENT,
    CANONICAL
}
