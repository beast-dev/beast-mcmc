package dr.inference.timeseries.engine.gmrf;

import dr.inference.model.GradientProvider;
import dr.inference.model.Parameter;
import dr.inference.timeseries.engine.GradientEngine;
import dr.math.distributions.GaussianMarkovRandomField;

/**
 * Gradient engine delegating to the existing GaussianMarkovRandomField implementation.
 *
 * It supports gradients with respect to the field itself and any parameter explicitly supported
 * by the wrapped GMRF via getGradientWrt(parameter).
 */
public class SparsePrecisionGradientEngine implements GradientEngine {

    private final GaussianMarkovRandomField gmrf;
    private final Parameter fieldParameter;

    public SparsePrecisionGradientEngine(final GaussianMarkovRandomField gmrf,
                                         final Parameter fieldParameter) {
        if (gmrf == null) {
            throw new IllegalArgumentException("gmrf must not be null");
        }
        if (fieldParameter == null) {
            throw new IllegalArgumentException("fieldParameter must not be null");
        }
        this.gmrf = gmrf;
        this.fieldParameter = fieldParameter;
    }

    @Override
    public boolean supportsGradientWrt(final Parameter parameter) {
        if (parameter == fieldParameter) {
            return true;
        }
        try {
            return gmrf.getGradientWrt(parameter) != null;
        } catch (RuntimeException ex) {
            return false;
        }
    }

    @Override
    public double[] getGradientWrt(final Parameter parameter) {
        if (parameter == fieldParameter) {
            return gmrf.getGradientLogDensity(fieldParameter.getParameterValues());
        }

        final GradientProvider provider = gmrf.getGradientWrt(parameter);
        return provider.getGradientLogDensity(parameter.getParameterValues());
    }

    @Override
    public void makeDirty() {
        // The wrapped GMRF manages its own cache invalidation.
    }
}
