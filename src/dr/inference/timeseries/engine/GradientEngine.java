package dr.inference.timeseries.engine;

import dr.inference.model.Parameter;

/**
 * Internal gradient source, intentionally independent of BEAST's HMC adapter interface.
 */
public interface GradientEngine {

    boolean supportsGradientWrt(Parameter parameter);

    double[] getGradientWrt(Parameter parameter);

    void makeDirty();
}
