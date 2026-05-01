package dr.evomodel.treedatalikelihood.continuous.canonical.message;

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

    // -----------------------------------------------------------------------
    // Builder
    // -----------------------------------------------------------------------

    /** Returns a {@link Builder} pre-populated with the values of this instance. */
    public Builder toBuilder() {
        return new Builder()
                .relativeJitter(relativeJitter)
                .absoluteJitter(absoluteJitter)
                .robustInversionAttempts(robustInversionAttempts)
                .strictInversionAttempts(strictInversionAttempts)
                .forceStrictSpdInversion(forceStrictSpdInversion)
                .spdFailureDumpEnabled(spdFailureDumpEnabled)
                .pivotFloorDebugEnabled(pivotFloorDebugEnabled);
    }

    /**
     * Fluent builder for {@link CanonicalNumericsOptions}.
     * Start from one of the predefined instances ({@link #OU_TREE}, {@link #ORTHOGONAL_BLOCK})
     * via {@link CanonicalNumericsOptions#toBuilder()}, or from defaults via {@code new Builder()}.
     */
    public static final class Builder {

        private double  relativeJitter             = 1.0e-12;
        private double  absoluteJitter             = 1.0e-12;
        private int     robustInversionAttempts    = 8;
        private int     strictInversionAttempts    = 12;
        private boolean forceStrictSpdInversion    = false;
        private boolean spdFailureDumpEnabled      = false;
        private boolean pivotFloorDebugEnabled     = false;

        public Builder() { }

        public Builder relativeJitter(double value)          { relativeJitter = value;          return this; }
        public Builder absoluteJitter(double value)          { absoluteJitter = value;          return this; }
        public Builder robustInversionAttempts(int value)    { robustInversionAttempts = value; return this; }
        public Builder strictInversionAttempts(int value)    { strictInversionAttempts = value; return this; }
        public Builder forceStrictSpdInversion(boolean value){ forceStrictSpdInversion = value; return this; }
        public Builder spdFailureDumpEnabled(boolean value)  { spdFailureDumpEnabled = value;   return this; }
        public Builder pivotFloorDebugEnabled(boolean value) { pivotFloorDebugEnabled = value;  return this; }

        public CanonicalNumericsOptions build() {
            return new CanonicalNumericsOptions(
                    relativeJitter,
                    absoluteJitter,
                    robustInversionAttempts,
                    strictInversionAttempts,
                    forceStrictSpdInversion,
                    spdFailureDumpEnabled,
                    pivotFloorDebugEnabled);
        }
    }
}
