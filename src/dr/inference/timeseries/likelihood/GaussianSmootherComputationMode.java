package dr.inference.timeseries.likelihood;

/**
 * Choice of smoother backend for Gaussian time-series inference.
 */
public enum GaussianSmootherComputationMode {
    EXPECTATION,
    CANONICAL
}
