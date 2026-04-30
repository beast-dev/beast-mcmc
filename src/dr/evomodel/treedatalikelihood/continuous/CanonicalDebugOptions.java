package dr.evomodel.treedatalikelihood.continuous;

/**
 * Debug switches for the canonical OU tree pathway.
 *
 * <p>The property names are intentionally kept compatible with the historical
 * ad hoc checks. Production code should depend on this typed snapshot instead
 * of reading system properties in hot paths.</p>
 */
public final class CanonicalDebugOptions {

    private static final String NONORTH_OMEGA_PROPERTY = "beast.debug.nonorth.omega";
    private static final String SELECTION_COMPONENTS_PROPERTY = "beast.debug.selectionComponents";
    private static final String BRANCH_LOCAL_SELECTION_FD_PROPERTY = "beast.debug.branchLocalSelectionFD";
    private static final String BRANCH_LOCAL_MEAN_FD_PROPERTY = "beast.debug.branchLocalMeanFD";
    private static final String BRANCH_LOCAL_MEAN_FD_NODE_PROPERTY = "beast.debug.branchLocalMeanFD.node";
    private static final String NATIVE_ASSEMBLY_PROPERTY = "beast.debug.branchGradient.nativeAssembly";
    private static final String NATIVE_ASSEMBLY_NODE_PROPERTY = "beast.debug.branchGradient.nativeAssembly.node";
    private static final String NONFINITE_BRANCH_STATS_PROPERTY = "beast.debug.ou.nonfiniteBranchStats";
    private static final String TRANSITION_CACHE_PROPERTY = "beast.debug.canonicalTransitionCache";
    private static final String ASSERT_NO_GRADIENT_CACHE_MISSES_PROPERTY =
            "beast.debug.canonicalAssertNoGradientCacheMisses";

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

    private CanonicalDebugOptions(final boolean nonOrthogonalOmegaEnabled,
                                  final boolean selectionComponentsEnabled,
                                  final boolean branchLocalSelectionFiniteDifferenceEnabled,
                                  final boolean branchLocalMeanFiniteDifferenceEnabled,
                                  final Integer branchLocalMeanFiniteDifferenceNode,
                                  final boolean nativeAssemblyEnabled,
                                  final Integer nativeAssemblyNode,
                                  final boolean nonFiniteBranchStatsEnabled,
                                  final boolean transitionCacheDiagnosticsEnabled,
                                  final boolean assertNoGradientCacheMissesEnabled) {
        this.nonOrthogonalOmegaEnabled = nonOrthogonalOmegaEnabled;
        this.selectionComponentsEnabled = selectionComponentsEnabled;
        this.branchLocalSelectionFiniteDifferenceEnabled = branchLocalSelectionFiniteDifferenceEnabled;
        this.branchLocalMeanFiniteDifferenceEnabled = branchLocalMeanFiniteDifferenceEnabled;
        this.branchLocalMeanFiniteDifferenceNode = branchLocalMeanFiniteDifferenceNode;
        this.nativeAssemblyEnabled = nativeAssemblyEnabled;
        this.nativeAssemblyNode = nativeAssemblyNode;
        this.nonFiniteBranchStatsEnabled = nonFiniteBranchStatsEnabled;
        this.transitionCacheDiagnosticsEnabled = transitionCacheDiagnosticsEnabled;
        this.assertNoGradientCacheMissesEnabled = assertNoGradientCacheMissesEnabled;
    }

    public static CanonicalDebugOptions fromSystemProperties() {
        return new CanonicalDebugOptions(
                Boolean.getBoolean(NONORTH_OMEGA_PROPERTY),
                Boolean.getBoolean(SELECTION_COMPONENTS_PROPERTY),
                Boolean.getBoolean(BRANCH_LOCAL_SELECTION_FD_PROPERTY),
                Boolean.getBoolean(BRANCH_LOCAL_MEAN_FD_PROPERTY),
                integerProperty(BRANCH_LOCAL_MEAN_FD_NODE_PROPERTY),
                Boolean.getBoolean(NATIVE_ASSEMBLY_PROPERTY),
                integerProperty(NATIVE_ASSEMBLY_NODE_PROPERTY),
                Boolean.getBoolean(NONFINITE_BRANCH_STATS_PROPERTY),
                Boolean.getBoolean(TRANSITION_CACHE_PROPERTY),
                Boolean.getBoolean(ASSERT_NO_GRADIENT_CACHE_MISSES_PROPERTY));
    }

    public boolean isNonOrthogonalOmegaEnabled() {
        return nonOrthogonalOmegaEnabled;
    }

    public boolean isSelectionComponentsEnabled() {
        return selectionComponentsEnabled;
    }

    public boolean isBranchLocalSelectionFiniteDifferenceEnabled() {
        return branchLocalSelectionFiniteDifferenceEnabled;
    }

    public boolean isBranchLocalMeanFiniteDifferenceEnabled() {
        return branchLocalMeanFiniteDifferenceEnabled;
    }

    public Integer getBranchLocalMeanFiniteDifferenceNode() {
        return branchLocalMeanFiniteDifferenceNode;
    }

    public boolean isNativeAssemblyEnabled() {
        return nativeAssemblyEnabled;
    }

    public Integer getNativeAssemblyNode() {
        return nativeAssemblyNode;
    }

    public boolean isNonFiniteBranchStatsEnabled() {
        return nonFiniteBranchStatsEnabled;
    }

    public boolean isTransitionCacheDiagnosticsEnabled() {
        return transitionCacheDiagnosticsEnabled;
    }

    public boolean isAssertNoGradientCacheMissesEnabled() {
        return assertNoGradientCacheMissesEnabled;
    }

    private static Integer integerProperty(final String property) {
        final String raw = System.getProperty(property);
        if (raw == null || raw.trim().isEmpty()) {
            return null;
        }
        return Integer.valueOf(raw.trim());
    }
}
