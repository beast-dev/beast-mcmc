package dr.evomodel.continuous.ou.canonical;

import dr.evomodel.continuous.ou.SelectionMatrixParameterization;

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
