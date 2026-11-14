package dr.evomodel.continuous.ou.math;

/**
 * Numerical policy for flat positive-definite matrix inversion.
 */
public interface PositiveDefiniteMatrixInversionOptions {

    double jitterBase(double maxAbsDiagonal);

    boolean isPivotFloorDebugEnabled();

    default void reportPivotFloor(final int clippedPivotCount,
                                  final int dimension,
                                  final double minPivotBeforeFloor,
                                  final double pivotFloor,
                                  final double maxAbsDiagonal) {
    }
}
