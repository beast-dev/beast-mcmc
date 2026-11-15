package dr.evomodel.continuous.ou.math;

/**
 * Standard OU-layer numerical policies for positive-definite matrix inversion.
 */
public final class OUPositiveDefiniteMatrixInversionOptions
        implements PositiveDefiniteMatrixInversionOptions {

    private static final String ORTHOGONAL_DEBUG_PD_PIVOT_FLOOR_PROPERTY =
            "beast.debug.pdPivotFloor";

    public static final OUPositiveDefiniteMatrixInversionOptions OU_TREE =
            new OUPositiveDefiniteMatrixInversionOptions(
                    1.0e-12,
                    1.0e-12,
                    false);

    public static final OUPositiveDefiniteMatrixInversionOptions ORTHOGONAL_BLOCK =
            new OUPositiveDefiniteMatrixInversionOptions(
                    1.0e-14,
                    1.0e-14,
                    Boolean.getBoolean(ORTHOGONAL_DEBUG_PD_PIVOT_FLOOR_PROPERTY));

    private final double relativeJitter;
    private final double absoluteJitter;
    private final boolean pivotFloorDebugEnabled;

    public OUPositiveDefiniteMatrixInversionOptions(final double relativeJitter,
                                                    final double absoluteJitter,
                                                    final boolean pivotFloorDebugEnabled) {
        this.relativeJitter = relativeJitter;
        this.absoluteJitter = absoluteJitter;
        this.pivotFloorDebugEnabled = pivotFloorDebugEnabled;
    }

    @Override
    public double jitterBase(final double maxAbsDiagonal) {
        return Math.max(absoluteJitter, relativeJitter * Math.max(1.0, maxAbsDiagonal));
    }

    @Override
    public boolean isPivotFloorDebugEnabled() {
        return pivotFloorDebugEnabled;
    }

    @Override
    public void reportPivotFloor(final int clippedPivotCount,
                                 final int dimension,
                                 final double minPivotBeforeFloor,
                                 final double pivotFloor,
                                 final double maxAbsDiagonal) {
        if (pivotFloorDebugEnabled && clippedPivotCount > 0) {
            System.err.println("pdPivotFloorDebug clipped=" + clippedPivotCount
                    + " dim=" + dimension
                    + " minPivotBeforeFloor=" + minPivotBeforeFloor
                    + " pivotFloor=" + pivotFloor
                    + " maxAbsDiagonal=" + maxAbsDiagonal);
        }
    }
}
