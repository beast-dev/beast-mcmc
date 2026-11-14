package dr.evomodel.treedatalikelihood.continuousNew;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

/**
 * Debug switches for the canonical OU tree pathway.
 *
 * <p>Create instances via {@link #fromSystemProperties()} (production) or
 * {@link #builder()} (tests and explicit configuration). The builder sets all
 * fields to the same defaults that {@code fromSystemProperties()} would use when
 * no system properties are set (i.e. everything disabled).
 *
 * <p>Call {@link #reportUnknownProperties()} at startup to surface typos in
 * property names before they silently have no effect.
 */
public final class CanonicalDebugOptions {

    // -----------------------------------------------------------------------
    // Known property names
    // -----------------------------------------------------------------------

    private static final String TRANSITION_CACHE_PROPERTY                 = "beast.debug.canonicalTransitionCache";
    private static final String ASSERT_NO_GRADIENT_CACHE_MISSES_PROPERTY  =
            "beast.debug.canonicalAssertNoGradientCacheMisses";

    /**
     * All {@code beast.debug.*} property names read anywhere in the canonical OU code,
     * including the ones read directly by the timers and numerics helpers rather than
     * through this class.  Used by {@link #reportUnknownProperties()} to detect typos.
     */
    public static final Set<String> KNOWN_PROPERTIES = Collections.unmodifiableSet(
            new HashSet<>(Arrays.asList(
                    TRANSITION_CACHE_PROPERTY,
                    "beast.debug.canonicalTransitionCacheReportEvery",
                    ASSERT_NO_GRADIENT_CACHE_MISSES_PROPERTY,
                    "beast.debug.canonicalPhaseTiming",
                    "beast.debug.canonicalPhaseTimingReportEvery",
                    "beast.debug.canonicalTraversalTiming",
                    "beast.debug.canonicalTraversalTimingReportEvery",
                    "beast.debug.canonicalFrechetPlan",
                    "beast.debug.canonicalFrechetPlanReportEvery",
                    "beast.debug.canonicalMatrixKernels",
                    "beast.debug.canonicalMatrixKernelsReportEvery",
                    "beast.debug.pdPivotFloor",
                    "beast.debug.ou.forceStrictSpdInverse",
                    "beast.debug.ou.spdDump")));

    // -----------------------------------------------------------------------
    // Fields
    // -----------------------------------------------------------------------

    private final boolean transitionCacheDiagnosticsEnabled;
    private final boolean assertNoGradientCacheMissesEnabled;

    private CanonicalDebugOptions(final Builder b) {
        this.transitionCacheDiagnosticsEnabled                 = b.transitionCacheDiagnosticsEnabled;
        this.assertNoGradientCacheMissesEnabled                = b.assertNoGradientCacheMissesEnabled;
    }

    // -----------------------------------------------------------------------
    // Factories
    // -----------------------------------------------------------------------

    /** Returns a {@link Builder} with all defaults (all flags off). */
    public static Builder builder() {
        return new Builder();
    }

    /**
     * Returns an instance populated from system properties.
     * Also calls {@link #reportUnknownProperties()} to warn about unrecognized
     * {@code beast.debug.*} property names.
     */
    public static CanonicalDebugOptions fromSystemProperties() {
        reportUnknownProperties();
        return builder()
                .transitionCacheDiagnosticsEnabled(Boolean.getBoolean(TRANSITION_CACHE_PROPERTY))
                .assertNoGradientCacheMissesEnabled(Boolean.getBoolean(ASSERT_NO_GRADIENT_CACHE_MISSES_PROPERTY))
                .build();
    }

    /**
     * Scans all current system properties for names that start with {@code "beast.debug."}
     * but are not in {@link #KNOWN_PROPERTIES}, and prints a warning for each.
     *
     * <p>Call this once at application startup or from {@link #fromSystemProperties()}.
     * Silent typos (e.g. {@code beast.debug.canonicalTransitionCachee}) will be reported here
     * rather than being silently ignored.
     */
    public static void reportUnknownProperties() {
        System.getProperties().stringPropertyNames().stream()
                .filter(name -> name.startsWith("beast.debug."))
                .filter(name -> !KNOWN_PROPERTIES.contains(name))
                .sorted()
                .forEach(name -> dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalDiagnosticsLog.warning(
                        "WARNING: Unrecognized canonical debug property '" + name
                        + "' — check spelling against CanonicalDebugOptions.KNOWN_PROPERTIES"));
    }

    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    public boolean isTransitionCacheDiagnosticsEnabled()             { return transitionCacheDiagnosticsEnabled; }
    public boolean isAssertNoGradientCacheMissesEnabled()            { return assertNoGradientCacheMissesEnabled; }

    // -----------------------------------------------------------------------
    // Builder
    // -----------------------------------------------------------------------

    /**
     * Fluent builder for {@link CanonicalDebugOptions}.
     * All fields default to disabled — identical to the behavior of
     * {@link #fromSystemProperties()} when no system properties are set.
     */
    public static final class Builder {

        private boolean transitionCacheDiagnosticsEnabled           = false;
        private boolean assertNoGradientCacheMissesEnabled          = false;

        private Builder() { }

        public Builder transitionCacheDiagnosticsEnabled(boolean v)           { transitionCacheDiagnosticsEnabled = v;           return this; }
        public Builder assertNoGradientCacheMissesEnabled(boolean v)          { assertNoGradientCacheMissesEnabled = v;          return this; }

        public CanonicalDebugOptions build() {
            return new CanonicalDebugOptions(this);
        }
    }
}
