package dr.inference.timeseries.representation;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.inference.timeseries.core.TimeGrid;

/**
 * Optional extension for transition representations that cache repeated time increments.
 */
public interface CachedGaussianTransitionRepresentation extends GaussianTransitionRepresentation {

    void prepareTimeGrid(TimeGrid timeGrid);

    void getCanonicalTransition(int fromIndex,
                                int toIndex,
                                TimeGrid timeGrid,
                                CanonicalGaussianTransition out);

    RepeatedDeltaCacheStatistics getCacheStatistics();

    void makeDirty();
}
