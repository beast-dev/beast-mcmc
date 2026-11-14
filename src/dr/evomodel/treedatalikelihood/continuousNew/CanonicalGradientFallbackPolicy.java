package dr.evomodel.treedatalikelihood.continuousNew;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

/**
 * Experimental and fallback switches for canonical OU gradients.
 *
 * <p>Create instances via {@link #fromSystemProperties()} (production) or
 * {@link #builder()} (tests and explicit configuration). The builder sets all
 * fields to the same defaults that {@code fromSystemProperties()} would use when
 * no system properties are set.
 *
 * <p>Call {@link #reportUnknownProperties()} at startup to surface typos in
 * property names before they silently have no effect.
 */
public final class CanonicalGradientFallbackPolicy {

    // -----------------------------------------------------------------------
    // Known property names
    // -----------------------------------------------------------------------

    private static final String BRANCH_LENGTH_FINITE_DIFFERENCE_PROPERTY  = "beast.experimental.canonicalBranchLengthGradientFiniteDifference";
    private static final String BRANCH_GRADIENT_THREADS_PROPERTY          = "beast.experimental.canonicalBranchGradientThreads";
    private static final String TRANSITION_PRELOAD_THREADS_PROPERTY        = "beast.experimental.canonicalTransitionPreloadThreads";
    private static final String BLOCK_BASIS_HOMOGENEOUS_PROPERTY          = CanonicalOUOptions.BLOCK_BASIS_HOMOGENEOUS_PROPERTY;

    /**
     * All property names recognized by this class, plus the canonical
     * {@code beast.experimental.*} properties that are read directly by other classes.
     * Used by {@link #reportUnknownProperties()} to detect typos.
     */
    public static final Set<String> KNOWN_PROPERTIES = Collections.unmodifiableSet(
            new HashSet<>(Arrays.asList(
                    BRANCH_LENGTH_FINITE_DIFFERENCE_PROPERTY,
                    BRANCH_GRADIENT_THREADS_PROPERTY,
                    TRANSITION_PRELOAD_THREADS_PROPERTY,
                    BLOCK_BASIS_HOMOGENEOUS_PROPERTY,
                    "beast.experimental.canonicalUseCholeskySolveMessageProduct")));

    // -----------------------------------------------------------------------
    // Fields
    // -----------------------------------------------------------------------

    private final boolean branchLengthFiniteDifference;
    private final int branchGradientParallelism;
    private final int transitionPreloadParallelism;

    private CanonicalGradientFallbackPolicy(final Builder b) {
        this.branchLengthFiniteDifference            = b.branchLengthFiniteDifference;
        this.branchGradientParallelism               = b.branchGradientParallelism;
        this.transitionPreloadParallelism            = b.transitionPreloadParallelism;
    }

    // -----------------------------------------------------------------------
    // Factories
    // -----------------------------------------------------------------------

    /** Returns a {@link Builder} with all defaults (all flags off, parallelism = 1). */
    public static Builder builder() {
        return new Builder();
    }

    /**
     * Returns an instance populated from system properties.
     * Also calls {@link #reportUnknownProperties()} to warn about unrecognized
     * {@code beast.experimental.*} property names.
     */
    public static CanonicalGradientFallbackPolicy fromSystemProperties() {
        reportUnknownProperties();
        final int branchGradientParallelism = readBranchGradientParallelism();
        return builder()
                .branchLengthFiniteDifference(Boolean.getBoolean(BRANCH_LENGTH_FINITE_DIFFERENCE_PROPERTY))
                .branchGradientParallelism(branchGradientParallelism)
                .transitionPreloadParallelism(readTransitionPreloadParallelism(branchGradientParallelism))
                .build();
    }

    /**
     * Scans all current system properties for names that start with
     * {@code "beast.experimental."} but are not in {@link #KNOWN_PROPERTIES},
     * and prints a warning for each.
     */
    public static void reportUnknownProperties() {
        System.getProperties().stringPropertyNames().stream()
                .filter(name -> name.startsWith("beast.experimental.canonical")
                             || name.startsWith("beast.experimental.localSelection")
                             || name.startsWith("beast.experimental.noTranspose")
                             || name.startsWith("beast.experimental.orthBlock")
                             || name.startsWith("beast.experimental.transposeNative")
                             || name.startsWith("beast.experimental.nativeFrom")
                             || name.startsWith("beast.experimental.disableExact"))
                .filter(name -> !KNOWN_PROPERTIES.contains(name))
                .sorted()
                .forEach(name -> dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalDiagnosticsLog.warning(
                        "WARNING: Unrecognized canonical gradient property '" + name
                        + "' — check spelling against CanonicalGradientFallbackPolicy.KNOWN_PROPERTIES"));
    }

    /**
     * Convenience accessor for the branch-gradient thread count read from system properties.
     *
     * @deprecated Prefer {@link #fromSystemProperties()} which populates
     *     {@link #getBranchGradientParallelism()} as part of the full config object.
     *     This static method is retained for call sites that only need the thread count.
     */
    @Deprecated
    public static int branchGradientParallelismFromSystemProperties() {
        return readBranchGradientParallelism();
    }

    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    public boolean useBranchLengthFiniteDifference()           { return branchLengthFiniteDifference; }

    /** Number of threads to use for parallel branch-gradient evaluation (>= 1). */
    public int getBranchGradientParallelism()                  { return branchGradientParallelism; }

    /** Number of threads to use for prewarming canonical branch transitions (>= 1). */
    public int getTransitionPreloadParallelism()                { return transitionPreloadParallelism; }

    // -----------------------------------------------------------------------
    // Builder
    // -----------------------------------------------------------------------

    /**
     * Fluent builder for {@link CanonicalGradientFallbackPolicy}.
     * All flags default to {@code false}; the parallelism settings default to 1.
     */
    public static final class Builder {

        private boolean branchLengthFiniteDifference             = false;
        private int     branchGradientParallelism                = 1;
        private int     transitionPreloadParallelism             = 1;

        private Builder() { }

        public Builder branchLengthFiniteDifference(boolean v)             { branchLengthFiniteDifference = v;             return this; }

        /** Thread count for parallel branch-gradient evaluation; must be >= 1. */
        public Builder branchGradientParallelism(int v) {
            if (v < 1) throw new IllegalArgumentException("branchGradientParallelism must be >= 1, got " + v);
            branchGradientParallelism = v;
            return this;
        }

        /** Thread count for transition preloading; must be >= 1. */
        public Builder transitionPreloadParallelism(int v) {
            if (v < 1) throw new IllegalArgumentException("transitionPreloadParallelism must be >= 1, got " + v);
            transitionPreloadParallelism = v;
            return this;
        }

        public CanonicalGradientFallbackPolicy build() {
            return new CanonicalGradientFallbackPolicy(this);
        }
    }

    // -----------------------------------------------------------------------
    // Private helpers
    // -----------------------------------------------------------------------

    private static int readBranchGradientParallelism() {
        final String value = System.getProperty(BRANCH_GRADIENT_THREADS_PROPERTY);
        if (value != null) {
            try {
                return Math.max(1, Integer.parseInt(value));
            } catch (NumberFormatException e) {
                dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalDiagnosticsLog.warning("WARNING: Canonical gradient property '"
                        + BRANCH_GRADIENT_THREADS_PROPERTY
                        + "' has non-integer value '" + value + "'; using 1.");
            }
        }
        return 1;
    }

    private static int readTransitionPreloadParallelism(final int branchGradientParallelism) {
        final String value = System.getProperty(TRANSITION_PRELOAD_THREADS_PROPERTY);
        if (value != null) {
            try {
                return Math.max(1, Integer.parseInt(value));
            } catch (NumberFormatException e) {
                dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalDiagnosticsLog.warning("WARNING: Canonical gradient property '"
                        + TRANSITION_PRELOAD_THREADS_PROPERTY
                        + "' has non-integer value '" + value + "'; using default.");
            }
        }
        return 1;
    }
}
