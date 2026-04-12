package dr.inference.timeseries.engine.gaussian;

/**
 * Choice of internal algebra for forward Gaussian likelihood evaluation.
 */
public enum GaussianForwardComputationMode {
    EXPECTATION,
    CANONICAL
}
