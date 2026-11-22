package dr.evomodel.treedatalikelihood.continuousNew.canonical.adapter;

import dr.evomodel.treedatalikelihood.continuousNew.CanonicalDebugOptions;

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
