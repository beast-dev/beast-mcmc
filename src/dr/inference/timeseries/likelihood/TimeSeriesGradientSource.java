package dr.inference.timeseries.likelihood;

import dr.inference.model.GradientProvider;
import dr.inference.model.Parameter;

/**
 * Likelihood-side access point for per-parameter gradient providers.
 */
public interface TimeSeriesGradientSource {

    GradientProvider getGradientWrt(Parameter parameter);
}
