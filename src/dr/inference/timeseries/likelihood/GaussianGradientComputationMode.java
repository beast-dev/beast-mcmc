package dr.inference.timeseries.likelihood;

/**
 * Choice of gradient backend for Gaussian time-series likelihoods.
 */
public enum GaussianGradientComputationMode {
    DISABLED,
    EXPECTATION_ANALYTICAL,
    CANONICAL_ANALYTICAL
}
