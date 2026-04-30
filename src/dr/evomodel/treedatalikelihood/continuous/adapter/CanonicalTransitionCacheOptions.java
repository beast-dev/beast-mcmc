package dr.evomodel.treedatalikelihood.continuous.adapter;

import dr.evomodel.treedatalikelihood.continuous.CanonicalDebugOptions;

final class CanonicalTransitionCacheOptions {

    private final boolean diagnosticsEnabled;

    private CanonicalTransitionCacheOptions(final boolean diagnosticsEnabled) {
        this.diagnosticsEnabled = diagnosticsEnabled;
    }

    static CanonicalTransitionCacheOptions fromSystemProperties() {
        final CanonicalDebugOptions debugOptions = CanonicalDebugOptions.fromSystemProperties();
        return new CanonicalTransitionCacheOptions(
                debugOptions.isTransitionCacheDiagnosticsEnabled()
                        || debugOptions.isAssertNoGradientCacheMissesEnabled());
    }

    boolean isDiagnosticsEnabled() {
        return diagnosticsEnabled;
    }
}
