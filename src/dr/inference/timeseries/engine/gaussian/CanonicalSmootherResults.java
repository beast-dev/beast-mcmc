package dr.inference.timeseries.engine.gaussian;

import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;

/**
 * Shared access surface for canonical smoother backends.
 */
public interface CanonicalSmootherResults {

    CanonicalForwardTrajectory getCanonicalTrajectory();

    CanonicalBranchGradientCache getCanonicalBranchGradientCache();

    GaussianTransitionRepresentation getTransitionRepresentation();

    TimeGrid getTimeGrid();

    void makeDirty();
}
