package dr.evomodel.treedatalikelihood.continuous;

/**
 * Experimental and fallback switches for canonical OU gradients.
 */
public final class CanonicalGradientFallbackPolicy {

    private static final String NUMERIC_LOCAL_SELECTION_PROPERTY =
            "beast.experimental.localSelectionNumeric";
    private static final String FROZEN_LOG_FACTOR_NUMERIC_PROPERTY =
            "beast.experimental.localSelectionNumeric.useFrozenLogFactor";
    private static final String NO_TRANSPOSE_DOMEGA_PROPERTY =
            "beast.experimental.noTransposeDOmega";
    private static final String DENSE_PULLBACK_FOR_ORTHOGONAL_PROPERTY =
            "beast.experimental.orthBlockUseDensePullback";
    private static final String TRANSPOSE_NATIVE_DF_PROPERTY =
            "beast.experimental.transposeNativeDF";
    private static final String NATIVE_FROM_CONTRIBUTION_PROPERTY =
            "beast.experimental.nativeFromContribution";
    private static final String BRANCH_LENGTH_FINITE_DIFFERENCE_PROPERTY =
            "beast.experimental.canonicalBranchLengthGradientFiniteDifference";
    private static final String DISABLE_EXACT_TIP_SHORTCUT_PROPERTY =
            "beast.experimental.disableExactTipShortcut";
    private static final String BRANCH_GRADIENT_THREADS_PROPERTY =
            "beast.experimental.canonicalBranchGradientThreads";

    private final boolean numericLocalSelectionGradient;
    private final boolean numericLocalSelectionFromFrozenLogFactor;
    private final boolean noTransposeDOmega;
    private final boolean densePullbackForOrthogonal;
    private final boolean transposeNativeDF;
    private final boolean nativeGradientFromContribution;
    private final boolean branchLengthFiniteDifference;
    private final boolean exactTipShortcutDisabled;

    private CanonicalGradientFallbackPolicy(final boolean numericLocalSelectionGradient,
                                            final boolean numericLocalSelectionFromFrozenLogFactor,
                                            final boolean noTransposeDOmega,
                                            final boolean densePullbackForOrthogonal,
                                            final boolean transposeNativeDF,
                                            final boolean nativeGradientFromContribution,
                                            final boolean branchLengthFiniteDifference,
                                            final boolean exactTipShortcutDisabled) {
        this.numericLocalSelectionGradient = numericLocalSelectionGradient;
        this.numericLocalSelectionFromFrozenLogFactor = numericLocalSelectionFromFrozenLogFactor;
        this.noTransposeDOmega = noTransposeDOmega;
        this.densePullbackForOrthogonal = densePullbackForOrthogonal;
        this.transposeNativeDF = transposeNativeDF;
        this.nativeGradientFromContribution = nativeGradientFromContribution;
        this.branchLengthFiniteDifference = branchLengthFiniteDifference;
        this.exactTipShortcutDisabled = exactTipShortcutDisabled;
    }

    public static CanonicalGradientFallbackPolicy fromSystemProperties() {
        return new CanonicalGradientFallbackPolicy(
                Boolean.getBoolean(NUMERIC_LOCAL_SELECTION_PROPERTY),
                Boolean.getBoolean(FROZEN_LOG_FACTOR_NUMERIC_PROPERTY),
                Boolean.getBoolean(NO_TRANSPOSE_DOMEGA_PROPERTY),
                Boolean.getBoolean(DENSE_PULLBACK_FOR_ORTHOGONAL_PROPERTY),
                Boolean.getBoolean(TRANSPOSE_NATIVE_DF_PROPERTY),
                Boolean.getBoolean(NATIVE_FROM_CONTRIBUTION_PROPERTY),
                Boolean.getBoolean(BRANCH_LENGTH_FINITE_DIFFERENCE_PROPERTY),
                Boolean.getBoolean(DISABLE_EXACT_TIP_SHORTCUT_PROPERTY));
    }

    public static int branchGradientParallelismFromSystemProperties() {
        final String propertyValue = System.getProperty(BRANCH_GRADIENT_THREADS_PROPERTY);
        if (propertyValue != null) {
            try {
                return Math.max(1, Integer.parseInt(propertyValue));
            } catch (NumberFormatException ignored) {
                // Fall back to the reproducible default below.
            }
        }
        return 1;
    }

    public boolean useNumericLocalSelectionGradient() {
        return numericLocalSelectionGradient;
    }

    public boolean useNumericLocalSelectionFromFrozenLogFactor() {
        return numericLocalSelectionFromFrozenLogFactor;
    }

    public boolean useNoTransposeDOmega() {
        return noTransposeDOmega;
    }

    public boolean useDensePullbackForOrthogonal() {
        return densePullbackForOrthogonal;
    }

    public boolean transposeNativeDF() {
        return transposeNativeDF;
    }

    public boolean useNativeGradientFromContribution() {
        return nativeGradientFromContribution;
    }

    public boolean useBranchLengthFiniteDifference() {
        return branchLengthFiniteDifference;
    }

    public boolean isExactTipShortcutDisabled() {
        return exactTipShortcutDisabled;
    }
}
