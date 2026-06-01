package dr.inference.timeseries.engine.kalman;

import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;

/**
 * Shared access surface for smoother backends that expose moment-form trajectories
 * and RTS statistics to analytical gradient formulas.
 */
public interface MomentSmootherResults {

    MomentBranchSmootherStats[] getSmootherStats();

    MomentForwardTrajectory getTrajectory();

    GaussianTransitionRepresentation getTransitionRepresentation();

    TimeGrid getTimeGrid();

    void makeDirty();
}
