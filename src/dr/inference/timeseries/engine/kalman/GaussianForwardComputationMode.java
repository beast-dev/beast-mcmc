package dr.inference.timeseries.engine.kalman;

/**
 * Choice of internal algebra for forward Gaussian likelihood evaluation.
 */
public enum GaussianForwardComputationMode {
    MOMENT,
    CANONICAL
}
