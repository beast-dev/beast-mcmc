package dr.evomodel.treedatalikelihood.continuous.canonical;

import dr.util.TaskPool;

/**
 * Optional capability for warming provider-owned canonical branch transitions.
 */
public interface CanonicalTransitionPreloader {

    void preloadCanonicalTransitions(int rootNodeIndex, TaskPool taskPool, int chunkSize);
}
