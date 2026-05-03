package dr.evomodel.treedatalikelihood.continuous;

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

    private static final String NONORTH_OMEGA_PROPERTY                    = "beast.debug.nonorth.omega";
    private static final String SELECTION_COMPONENTS_PROPERTY              = "beast.debug.selectionComponents";
    private static final String BRANCH_LOCAL_SELECTION_FD_PROPERTY        = "beast.debug.branchLocalSelectionFD";
    private static final String BRANCH_LOCAL_MEAN_FD_PROPERTY             = "beast.debug.branchLocalMeanFD";
    private static final String BRANCH_LOCAL_MEAN_FD_NODE_PROPERTY        = "beast.debug.branchLocalMeanFD.node";
    private static final String NATIVE_ASSEMBLY_PROPERTY                   = "beast.debug.branchGradient.nativeAssembly";
    private static final String NATIVE_ASSEMBLY_NODE_PROPERTY              = "beast.debug.branchGradient.nativeAssembly.node";
    private static final String NONFINITE_BRANCH_STATS_PROPERTY           = "beast.debug.ou.nonfiniteBranchStats";
    private static final String TRANSITION_CACHE_PROPERTY                 = "beast.debug.canonicalTransitionCache";
    private static final String ASSERT_NO_GRADIENT_CACHE_MISSES_PROPERTY  =
            "beast.debug.canonicalAssertNoGradientCacheMisses";
    private static final String CANONICAL_PHASE_TIMING_PROPERTY           =
            "beast.debug.canonicalPhaseTiming";
    private static final String CANONICAL_PHASE_TIMING_REPORT_EVERY_PROPERTY =
            "beast.debug.canonicalPhaseTimingReportEvery";
    private static final String CANONICAL_TRAVERSAL_TIMING_PROPERTY       =
            "beast.debug.canonicalTraversalTiming";
    private static final String CANONICAL_TRAVERSAL_TIMING_REPORT_EVERY_PROPERTY =
            "beast.debug.canonicalTraversalTimingReportEvery";
    private static final String CANONICAL_FRECHET_PLAN_PROPERTY           =
            "beast.debug.canonicalFrechetPlan";
    private static final String CANONICAL_FRECHET_PLAN_REPORT_EVERY_PROPERTY =
            "beast.debug.canonicalFrechetPlanReportEvery";

    /**
     * All property names recognized by this class.
     * Used by {@link #reportUnknownProperties()} to detect typos.
     */
    public static final Set<String> KNOWN_PROPERTIES = Collections.unmodifiableSet(
            new HashSet<>(Arrays.asList(
                    NONORTH_OMEGA_PROPERTY,
                    SELECTION_COMPONENTS_PROPERTY,
                    BRANCH_LOCAL_SELECTION_FD_PROPERTY,
                    BRANCH_LOCAL_MEAN_FD_PROPERTY,
                    BRANCH_LOCAL_MEAN_FD_NODE_PROPERTY,
                    NATIVE_ASSEMBLY_PROPERTY,
                    NATIVE_ASSEMBLY_NODE_PROPERTY,
                    NONFINITE_BRANCH_STATS_PROPERTY,
                    TRANSITION_CACHE_PROPERTY,
                    ASSERT_NO_GRADIENT_CACHE_MISSES_PROPERTY,
                    CANONICAL_PHASE_TIMING_PROPERTY,
                    CANONICAL_PHASE_TIMING_REPORT_EVERY_PROPERTY,
                    CANONICAL_TRAVERSAL_TIMING_PROPERTY,
                    CANONICAL_TRAVERSAL_TIMING_REPORT_EVERY_PROPERTY,
                    CANONICAL_FRECHET_PLAN_PROPERTY,
                    CANONICAL_FRECHET_PLAN_REPORT_EVERY_PROPERTY)));

    // -----------------------------------------------------------------------
    // Fields
    // -----------------------------------------------------------------------

    private final boolean nonOrthogonalOmegaEnabled;
    private final boolean selectionComponentsEnabled;
    private final boolean branchLocalSelectionFiniteDifferenceEnabled;
    private final boolean branchLocalMeanFiniteDifferenceEnabled;
    private final Integer branchLocalMeanFiniteDifferenceNode;
    private final boolean nativeAssemblyEnabled;
    private final Integer nativeAssemblyNode;
    private final boolean nonFiniteBranchStatsEnabled;
    private final boolean transitionCacheDiagnosticsEnabled;
    private final boolean assertNoGradientCacheMissesEnabled;

    private CanonicalDebugOptions(final Builder b) {
        this.nonOrthogonalOmegaEnabled                        = b.nonOrthogonalOmegaEnabled;
        this.selectionComponentsEnabled                        = b.selectionComponentsEnabled;
        this.branchLocalSelectionFiniteDifferenceEnabled       = b.branchLocalSelectionFiniteDifferenceEnabled;
        this.branchLocalMeanFiniteDifferenceEnabled            = b.branchLocalMeanFiniteDifferenceEnabled;
        this.branchLocalMeanFiniteDifferenceNode               = b.branchLocalMeanFiniteDifferenceNode;
        this.nativeAssemblyEnabled                             = b.nativeAssemblyEnabled;
        this.nativeAssemblyNode                                = b.nativeAssemblyNode;
        this.nonFiniteBranchStatsEnabled                       = b.nonFiniteBranchStatsEnabled;
        this.transitionCacheDiagnosticsEnabled                 = b.transitionCacheDiagnosticsEnabled;
        this.assertNoGradientCacheMissesEnabled                = b.assertNoGradientCacheMissesEnabled;
    }

    // -----------------------------------------------------------------------
    // Factories
    // -----------------------------------------------------------------------

    /** Returns a {@link Builder} with all defaults (all flags off, all nodes null). */
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
                .nonOrthogonalOmegaEnabled(Boolean.getBoolean(NONORTH_OMEGA_PROPERTY))
                .selectionComponentsEnabled(Boolean.getBoolean(SELECTION_COMPONENTS_PROPERTY))
                .branchLocalSelectionFiniteDifferenceEnabled(Boolean.getBoolean(BRANCH_LOCAL_SELECTION_FD_PROPERTY))
                .branchLocalMeanFiniteDifferenceEnabled(Boolean.getBoolean(BRANCH_LOCAL_MEAN_FD_PROPERTY))
                .branchLocalMeanFiniteDifferenceNode(integerProperty(BRANCH_LOCAL_MEAN_FD_NODE_PROPERTY))
                .nativeAssemblyEnabled(Boolean.getBoolean(NATIVE_ASSEMBLY_PROPERTY))
                .nativeAssemblyNode(integerProperty(NATIVE_ASSEMBLY_NODE_PROPERTY))
                .nonFiniteBranchStatsEnabled(Boolean.getBoolean(NONFINITE_BRANCH_STATS_PROPERTY))
                .transitionCacheDiagnosticsEnabled(Boolean.getBoolean(TRANSITION_CACHE_PROPERTY))
                .assertNoGradientCacheMissesEnabled(Boolean.getBoolean(ASSERT_NO_GRADIENT_CACHE_MISSES_PROPERTY))
                .build();
    }

    /**
     * Scans all current system properties for names that start with {@code "beast.debug."}
     * but are not in {@link #KNOWN_PROPERTIES}, and prints a warning for each.
     *
     * <p>Call this once at application startup or from {@link #fromSystemProperties()}.
     * Silent typos (e.g. {@code beast.debug.nonorth.omegaa}) will be reported here
     * rather than being silently ignored.
     */
    public static void reportUnknownProperties() {
        System.getProperties().stringPropertyNames().stream()
                .filter(name -> name.startsWith("beast.debug."))
                .filter(name -> !KNOWN_PROPERTIES.contains(name))
                .sorted()
                .forEach(name -> dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalDiagnosticsLog.warning(
                        "WARNING: Unrecognized canonical debug property '" + name
                        + "' — check spelling against CanonicalDebugOptions.KNOWN_PROPERTIES"));
    }

    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    public boolean isNonOrthogonalOmegaEnabled()                    { return nonOrthogonalOmegaEnabled; }
    public boolean isSelectionComponentsEnabled()                    { return selectionComponentsEnabled; }
    public boolean isBranchLocalSelectionFiniteDifferenceEnabled()   { return branchLocalSelectionFiniteDifferenceEnabled; }
    public boolean isBranchLocalMeanFiniteDifferenceEnabled()        { return branchLocalMeanFiniteDifferenceEnabled; }
    public Integer getBranchLocalMeanFiniteDifferenceNode()          { return branchLocalMeanFiniteDifferenceNode; }
    public boolean isNativeAssemblyEnabled()                         { return nativeAssemblyEnabled; }
    public Integer getNativeAssemblyNode()                           { return nativeAssemblyNode; }
    public boolean isNonFiniteBranchStatsEnabled()                   { return nonFiniteBranchStatsEnabled; }
    public boolean isTransitionCacheDiagnosticsEnabled()             { return transitionCacheDiagnosticsEnabled; }
    public boolean isAssertNoGradientCacheMissesEnabled()            { return assertNoGradientCacheMissesEnabled; }

    // -----------------------------------------------------------------------
    // Builder
    // -----------------------------------------------------------------------

    /**
     * Fluent builder for {@link CanonicalDebugOptions}.
     * All fields default to disabled / null — identical to the behavior of
     * {@link #fromSystemProperties()} when no system properties are set.
     */
    public static final class Builder {

        private boolean nonOrthogonalOmegaEnabled                  = false;
        private boolean selectionComponentsEnabled                  = false;
        private boolean branchLocalSelectionFiniteDifferenceEnabled = false;
        private boolean branchLocalMeanFiniteDifferenceEnabled      = false;
        private Integer branchLocalMeanFiniteDifferenceNode         = null;
        private boolean nativeAssemblyEnabled                       = false;
        private Integer nativeAssemblyNode                          = null;
        private boolean nonFiniteBranchStatsEnabled                 = false;
        private boolean transitionCacheDiagnosticsEnabled           = false;
        private boolean assertNoGradientCacheMissesEnabled          = false;

        private Builder() { }

        public Builder nonOrthogonalOmegaEnabled(boolean v)                  { nonOrthogonalOmegaEnabled = v;                  return this; }
        public Builder selectionComponentsEnabled(boolean v)                  { selectionComponentsEnabled = v;                  return this; }
        public Builder branchLocalSelectionFiniteDifferenceEnabled(boolean v) { branchLocalSelectionFiniteDifferenceEnabled = v; return this; }
        public Builder branchLocalMeanFiniteDifferenceEnabled(boolean v)      { branchLocalMeanFiniteDifferenceEnabled = v;      return this; }
        public Builder branchLocalMeanFiniteDifferenceNode(Integer v)         { branchLocalMeanFiniteDifferenceNode = v;         return this; }
        public Builder nativeAssemblyEnabled(boolean v)                       { nativeAssemblyEnabled = v;                       return this; }
        public Builder nativeAssemblyNode(Integer v)                          { nativeAssemblyNode = v;                          return this; }
        public Builder nonFiniteBranchStatsEnabled(boolean v)                 { nonFiniteBranchStatsEnabled = v;                 return this; }
        public Builder transitionCacheDiagnosticsEnabled(boolean v)           { transitionCacheDiagnosticsEnabled = v;           return this; }
        public Builder assertNoGradientCacheMissesEnabled(boolean v)          { assertNoGradientCacheMissesEnabled = v;          return this; }

        public CanonicalDebugOptions build() {
            return new CanonicalDebugOptions(this);
        }
    }

    // -----------------------------------------------------------------------
    // Private helpers
    // -----------------------------------------------------------------------

    private static Integer integerProperty(final String property) {
        final String raw = System.getProperty(property);
        if (raw == null || raw.trim().isEmpty()) {
            return null;
        }
        try {
            return Integer.valueOf(raw.trim());
        } catch (NumberFormatException e) {
            dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalDiagnosticsLog.warning("WARNING: Canonical debug property '" + property
                    + "' has non-integer value '" + raw + "'; ignoring.");
            return null;
        }
    }
}
