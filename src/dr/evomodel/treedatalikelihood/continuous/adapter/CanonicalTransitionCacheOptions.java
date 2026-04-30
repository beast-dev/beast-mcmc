package dr.evomodel.treedatalikelihood.continuous.adapter;

final class CanonicalTransitionCacheOptions {

    private static final String DEBUG_CACHE_PROPERTY = "beast.debug.canonicalTransitionCache";
    private static final String ASSERT_NO_GRADIENT_CACHE_MISSES_PROPERTY =
            "beast.debug.canonicalAssertNoGradientCacheMisses";

    private final boolean diagnosticsEnabled;

    private CanonicalTransitionCacheOptions(final boolean diagnosticsEnabled) {
        this.diagnosticsEnabled = diagnosticsEnabled;
    }

    static CanonicalTransitionCacheOptions fromSystemProperties() {
        return new CanonicalTransitionCacheOptions(
                Boolean.getBoolean(DEBUG_CACHE_PROPERTY)
                        || Boolean.getBoolean(ASSERT_NO_GRADIENT_CACHE_MISSES_PROPERTY));
    }

    boolean isDiagnosticsEnabled() {
        return diagnosticsEnabled;
    }
}
