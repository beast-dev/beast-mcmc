package dr.evomodel.treedatalikelihood.continuous.gaussian.message;

/**
 * Numerical policy knobs for canonical Gaussian/OU code.
 */
public final class CanonicalNumericsOptions {

    private static final String OU_SPD_DEBUG_DUMP_PROPERTY =
            "beast.debug.ou.spdDump";
    private static final String OU_FORCE_STRICT_SPD_INVERSION_PROPERTY =
            "beast.debug.ou.forceStrictSpdInverse";
    private static final String ORTHOGONAL_DEBUG_PD_PIVOT_FLOOR_PROPERTY =
            "beast.debug.pdPivotFloor";

    public static final CanonicalNumericsOptions OU_TREE = new CanonicalNumericsOptions(
            1.0e-12,
            1.0e-12,
            8,
            12,
            Boolean.getBoolean(OU_FORCE_STRICT_SPD_INVERSION_PROPERTY),
            Boolean.getBoolean(OU_SPD_DEBUG_DUMP_PROPERTY),
            false);

    public static final CanonicalNumericsOptions ORTHOGONAL_BLOCK = new CanonicalNumericsOptions(
            1.0e-14,
            1.0e-14,
            8,
            12,
            false,
            false,
            Boolean.getBoolean(ORTHOGONAL_DEBUG_PD_PIVOT_FLOOR_PROPERTY));

    private final double relativeJitter;
    private final double absoluteJitter;
    private final int robustInversionAttempts;
    private final int strictInversionAttempts;
    private final boolean forceStrictSpdInversion;
    private final boolean spdFailureDumpEnabled;
    private final boolean pivotFloorDebugEnabled;

    public CanonicalNumericsOptions(final double relativeJitter,
                                    final double absoluteJitter,
                                    final int robustInversionAttempts,
                                    final int strictInversionAttempts,
                                    final boolean forceStrictSpdInversion,
                                    final boolean spdFailureDumpEnabled,
                                    final boolean pivotFloorDebugEnabled) {
        this.relativeJitter = relativeJitter;
        this.absoluteJitter = absoluteJitter;
        this.robustInversionAttempts = robustInversionAttempts;
        this.strictInversionAttempts = strictInversionAttempts;
        this.forceStrictSpdInversion = forceStrictSpdInversion;
        this.spdFailureDumpEnabled = spdFailureDumpEnabled;
        this.pivotFloorDebugEnabled = pivotFloorDebugEnabled;
    }

    public double getRelativeJitter() {
        return relativeJitter;
    }

    public double getAbsoluteJitter() {
        return absoluteJitter;
    }

    public int getRobustInversionAttempts() {
        return robustInversionAttempts;
    }

    public int getStrictInversionAttempts() {
        return strictInversionAttempts;
    }

    public boolean isForceStrictSpdInversion() {
        return forceStrictSpdInversion;
    }

    public boolean isSpdFailureDumpEnabled() {
        return spdFailureDumpEnabled;
    }

    public boolean isPivotFloorDebugEnabled() {
        return pivotFloorDebugEnabled;
    }

    public double jitterBase(final double maxAbsDiagonal) {
        return Math.max(absoluteJitter, relativeJitter * Math.max(1.0, maxAbsDiagonal));
    }
}
