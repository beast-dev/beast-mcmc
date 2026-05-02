package dr.evomodel.treedatalikelihood.continuous;

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

    private static final String NUMERIC_LOCAL_SELECTION_PROPERTY          = "beast.experimental.localSelectionNumeric";
    private static final String FROZEN_LOG_FACTOR_NUMERIC_PROPERTY        = "beast.experimental.localSelectionNumeric.useFrozenLogFactor";
    private static final String NO_TRANSPOSE_DOMEGA_PROPERTY              = "beast.experimental.noTransposeDOmega";
    private static final String DENSE_PULLBACK_FOR_ORTHOGONAL_PROPERTY    = "beast.experimental.orthBlockUseDensePullback";
    private static final String TRANSPOSE_NATIVE_DF_PROPERTY              = "beast.experimental.transposeNativeDF";
    private static final String NATIVE_FROM_CONTRIBUTION_PROPERTY         = "beast.experimental.nativeFromContribution";
    private static final String BRANCH_LENGTH_FINITE_DIFFERENCE_PROPERTY  = "beast.experimental.canonicalBranchLengthGradientFiniteDifference";
    private static final String DISABLE_EXACT_TIP_SHORTCUT_PROPERTY       = "beast.experimental.disableExactTipShortcut";
    private static final String BRANCH_GRADIENT_THREADS_PROPERTY          = "beast.experimental.canonicalBranchGradientThreads";
    private static final String PREPARED_TRANSITION_MOMENTS_PROPERTY      = "beast.experimental.canonicalUsePreparedTransitionMoments";

    /**
     * All property names recognized by this class.
     * Used by {@link #reportUnknownProperties()} to detect typos.
     */
    public static final Set<String> KNOWN_PROPERTIES = Collections.unmodifiableSet(
            new HashSet<>(Arrays.asList(
                    NUMERIC_LOCAL_SELECTION_PROPERTY,
                    FROZEN_LOG_FACTOR_NUMERIC_PROPERTY,
                    NO_TRANSPOSE_DOMEGA_PROPERTY,
                    DENSE_PULLBACK_FOR_ORTHOGONAL_PROPERTY,
                    TRANSPOSE_NATIVE_DF_PROPERTY,
                    NATIVE_FROM_CONTRIBUTION_PROPERTY,
                    BRANCH_LENGTH_FINITE_DIFFERENCE_PROPERTY,
                    DISABLE_EXACT_TIP_SHORTCUT_PROPERTY,
                    BRANCH_GRADIENT_THREADS_PROPERTY,
                    PREPARED_TRANSITION_MOMENTS_PROPERTY)));

    // -----------------------------------------------------------------------
    // Fields
    // -----------------------------------------------------------------------

    private final boolean numericLocalSelectionGradient;
    private final boolean numericLocalSelectionFromFrozenLogFactor;
    private final boolean noTransposeDOmega;
    private final boolean densePullbackForOrthogonal;
    private final boolean transposeNativeDF;
    private final boolean nativeGradientFromContribution;
    private final boolean branchLengthFiniteDifference;
    private final boolean exactTipShortcutDisabled;
    private final int branchGradientParallelism;

    private CanonicalGradientFallbackPolicy(final Builder b) {
        this.numericLocalSelectionGradient           = b.numericLocalSelectionGradient;
        this.numericLocalSelectionFromFrozenLogFactor= b.numericLocalSelectionFromFrozenLogFactor;
        this.noTransposeDOmega                       = b.noTransposeDOmega;
        this.densePullbackForOrthogonal              = b.densePullbackForOrthogonal;
        this.transposeNativeDF                       = b.transposeNativeDF;
        this.nativeGradientFromContribution          = b.nativeGradientFromContribution;
        this.branchLengthFiniteDifference            = b.branchLengthFiniteDifference;
        this.exactTipShortcutDisabled                = b.exactTipShortcutDisabled;
        this.branchGradientParallelism               = b.branchGradientParallelism;
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
        return builder()
                .numericLocalSelectionGradient(Boolean.getBoolean(NUMERIC_LOCAL_SELECTION_PROPERTY))
                .numericLocalSelectionFromFrozenLogFactor(Boolean.getBoolean(FROZEN_LOG_FACTOR_NUMERIC_PROPERTY))
                .noTransposeDOmega(Boolean.getBoolean(NO_TRANSPOSE_DOMEGA_PROPERTY))
                .densePullbackForOrthogonal(Boolean.getBoolean(DENSE_PULLBACK_FOR_ORTHOGONAL_PROPERTY))
                .transposeNativeDF(Boolean.getBoolean(TRANSPOSE_NATIVE_DF_PROPERTY))
                .nativeGradientFromContribution(Boolean.getBoolean(NATIVE_FROM_CONTRIBUTION_PROPERTY))
                .branchLengthFiniteDifference(Boolean.getBoolean(BRANCH_LENGTH_FINITE_DIFFERENCE_PROPERTY))
                .exactTipShortcutDisabled(Boolean.getBoolean(DISABLE_EXACT_TIP_SHORTCUT_PROPERTY))
                .branchGradientParallelism(readParallelism())
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
                .forEach(name -> dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalDiagnosticsLog.warning(
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
        return readParallelism();
    }

    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    public boolean useNumericLocalSelectionGradient()          { return numericLocalSelectionGradient; }
    public boolean useNumericLocalSelectionFromFrozenLogFactor(){ return numericLocalSelectionFromFrozenLogFactor; }
    public boolean useNoTransposeDOmega()                      { return noTransposeDOmega; }
    public boolean useDensePullbackForOrthogonal()             { return densePullbackForOrthogonal; }
    public boolean transposeNativeDF()                         { return transposeNativeDF; }
    public boolean useNativeGradientFromContribution()         { return nativeGradientFromContribution; }
    public boolean useBranchLengthFiniteDifference()           { return branchLengthFiniteDifference; }
    public boolean isExactTipShortcutDisabled()                { return exactTipShortcutDisabled; }

    /** Number of threads to use for parallel branch-gradient evaluation (>= 1). */
    public int getBranchGradientParallelism()                  { return branchGradientParallelism; }

    // -----------------------------------------------------------------------
    // Builder
    // -----------------------------------------------------------------------

    /**
     * Fluent builder for {@link CanonicalGradientFallbackPolicy}.
     * All flags default to {@code false}; {@code branchGradientParallelism} defaults to 1.
     */
    public static final class Builder {

        private boolean numericLocalSelectionGradient            = false;
        private boolean numericLocalSelectionFromFrozenLogFactor = false;
        private boolean noTransposeDOmega                        = false;
        private boolean densePullbackForOrthogonal               = false;
        private boolean transposeNativeDF                        = false;
        private boolean nativeGradientFromContribution           = false;
        private boolean branchLengthFiniteDifference             = false;
        private boolean exactTipShortcutDisabled                 = false;
        private int     branchGradientParallelism                = 1;

        private Builder() { }

        public Builder numericLocalSelectionGradient(boolean v)            { numericLocalSelectionGradient = v;            return this; }
        public Builder numericLocalSelectionFromFrozenLogFactor(boolean v) { numericLocalSelectionFromFrozenLogFactor = v; return this; }
        public Builder noTransposeDOmega(boolean v)                        { noTransposeDOmega = v;                        return this; }
        public Builder densePullbackForOrthogonal(boolean v)               { densePullbackForOrthogonal = v;               return this; }
        public Builder transposeNativeDF(boolean v)                        { transposeNativeDF = v;                        return this; }
        public Builder nativeGradientFromContribution(boolean v)           { nativeGradientFromContribution = v;           return this; }
        public Builder branchLengthFiniteDifference(boolean v)             { branchLengthFiniteDifference = v;             return this; }
        public Builder exactTipShortcutDisabled(boolean v)                 { exactTipShortcutDisabled = v;                 return this; }

        /** Thread count for parallel branch-gradient evaluation; must be >= 1. */
        public Builder branchGradientParallelism(int v) {
            if (v < 1) throw new IllegalArgumentException("branchGradientParallelism must be >= 1, got " + v);
            branchGradientParallelism = v;
            return this;
        }

        public CanonicalGradientFallbackPolicy build() {
            return new CanonicalGradientFallbackPolicy(this);
        }
    }

    // -----------------------------------------------------------------------
    // Private helpers
    // -----------------------------------------------------------------------

    private static int readParallelism() {
        final String value = System.getProperty(BRANCH_GRADIENT_THREADS_PROPERTY);
        if (value != null) {
            try {
                return Math.max(1, Integer.parseInt(value));
            } catch (NumberFormatException e) {
                dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalDiagnosticsLog.warning("WARNING: Canonical gradient property '"
                        + BRANCH_GRADIENT_THREADS_PROPERTY
                        + "' has non-integer value '" + value + "'; using 1.");
            }
        }
        return 1;
    }
}
