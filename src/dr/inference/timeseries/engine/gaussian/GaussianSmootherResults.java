package dr.inference.timeseries.engine.gaussian;

import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;

/**
 * Shared access surface for smoother backends that expose moment-form trajectories
 * and RTS statistics to analytical gradient formulas.
 */
public interface GaussianSmootherResults {

    BranchSmootherStats[] getSmootherStats();

    ForwardTrajectory getTrajectory();

    GaussianTransitionRepresentation getTransitionRepresentation();

    TimeGrid getTimeGrid();

    void makeDirty();
}
