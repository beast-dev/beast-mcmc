package dr.inference.model;

import dr.inference.model.Likelihood;

/**
 * An interface for a prefetch-compatible likelihood.
 *
 * The calling operator calls startPrefetchOperation to specify which prefetch buffer is being targeted.
 * It then operates on the parameter, the likelihood is listening and records the change. Then
 * prefetchLogLikelihoods is called and all the updated likelihoods are calculated in parallel.
 * Finally acceptPrefetch is called to specify which (if any) of the prefetches has been accepted -
 * this may need to update the other buffers to reflect this.
 *
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface PrefetchableLikelihood extends Likelihood {
    /**
     * @return the number of prefetch operations that will be done in parallel
     */
    int getPrefetchCount();

    /**
     * This is called before the particular prefetch operation is done so the likelihood can handle
     * the operation appropriately (i.e., listening for changes).
     * @param prefetch
     */
    void startPrefetchOperation(int prefetch);

    /**
     * Called after a particular prefetch operation is done so the likelihood can clean up and store
     * necessary info.
     * @param prefetch
     */
    void finishPrefetchOperation(int prefetch);

    /**
     * This is called to specify which cached prefetch likelihood will be return at the next
     * getLogLikelihood call.
     * @param prefetch
     */
    void setPrefetchLikelihood(int prefetch);

    /**
     * Calculate the likelihood for all the prefetch operations in parallel
     */
    void prefetchLogLikelihoods();

    /**
     * Specifies which of the prefetch operations was successful (if any). This allows the likelihood
     * to update its current state to the successful one.
     * @param prefetch
     */
    void acceptPrefetch(int prefetch);

    /**
     * Called if all the prefetch operations were rejected. The likelihood should go back to its initial
     * state. Note that the standard store/restore methods will be called as each prefetched operation
     * is evaluated. It is likely that these will be ignored until this method is called.
     */
    void rejectAllPrefetches();

    /**
     * Turns off prefetch likelihood calculations (i.e., for a non-prefetch operator or a full evaluation).
     */
    void suspendPrefetch();

    /**
     * Release the used buffer indices - used for debugging purposes when doing a non-parallel prefetch.
     * @param prefetch
     */
    void releaseBufferIndices(int prefetch);

    void setIgnoreTreeEvents(boolean ignoreTreeEvents);
}
