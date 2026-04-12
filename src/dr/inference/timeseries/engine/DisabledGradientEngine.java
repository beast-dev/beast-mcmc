package dr.inference.timeseries.engine;

import dr.inference.model.Parameter;

/**
 * Gradient engine placeholder used when a likelihood is intentionally built without
 * gradient support.
 */
public final class DisabledGradientEngine implements GradientEngine {

    @Override
    public boolean supportsGradientWrt(final Parameter parameter) {
        return false;
    }

    @Override
    public double[] getGradientWrt(final Parameter parameter) {
        throw new IllegalArgumentException("Gradient support is disabled for this time-series likelihood");
    }

    @Override
    public void makeDirty() {
        // no-op
    }
}
