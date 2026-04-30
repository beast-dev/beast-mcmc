package dr.evomodel.continuous.ou;

/**
 * Optional canonical OU backend for selection parametrizations with a faster
 * native transition, covariance, or gradient path.
 */
public interface SpecializedCanonicalSelectionParameterization extends SelectionMatrixParameterization,
        CanonicalGradientPackingCapability,
        CanonicalPreparedTransitionCapability,
        CanonicalNativeSelectionGradientCapability,
        CanonicalNativeMeanGradientCapability,
        CanonicalNativeDiffusionGradientCapability {
}
