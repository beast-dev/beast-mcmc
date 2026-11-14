package dr.evomodel.treedatalikelihood.continuousNew.canonical;

import dr.util.TaskPool;

/**
 * Optional capability for warming provider-owned canonical branch transitions.
 */
public interface CanonicalTransitionPreloader {

    void preloadCanonicalTransitions(int rootNodeIndex, TaskPool taskPool);
}
