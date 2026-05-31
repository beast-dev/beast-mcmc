package dr.inference.timeseries.engine.kalman.formula;

import dr.inference.model.Parameter;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.engine.kalman.CanonicalBranchGradientCache;
import dr.inference.timeseries.engine.kalman.CanonicalForwardTrajectory;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;

/**
 * Canonical-path analogue of {@link GradientFormula}.
 */
public interface CanonicalGradientFormula {

    boolean supportsParameter(Parameter parameter);

    double[] computeGradient(Parameter parameter,
                             CanonicalForwardTrajectory trajectory,
                             GaussianTransitionRepresentation repr,
                             TimeGrid timeGrid);

    default double[] computeGradient(Parameter parameter,
                                     CanonicalForwardTrajectory trajectory,
                                     CanonicalBranchGradientCache branchGradientCache,
                                     GaussianTransitionRepresentation repr,
                                     TimeGrid timeGrid) {
        return computeGradient(parameter, trajectory, repr, timeGrid);
    }

    default void makeDirty() {
        // no-op
    }
}
