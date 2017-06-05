package dr.inference.model;

import dr.inference.model.Likelihood;

/**
 * An interface for a prefetch-compatible likelihood.
 *
 * The calling operator calls setCurrentPrefetch to specify which prefetch buffer is being targeted.
 * It then operates on the parameter, the likelihood is listening and records the change. Then
 * prefetchLogLikelihoods is called and all the updated likelihoods are calculated in parallel.
 * Finally acceptPrefetch is called to specify which (if any) of the prefetches has been accepted -
 * this may need to update the other buffers to reflect this.
 *
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface PrefetchableLikelihood extends Likelihood {
    int getPrefetchCount();

    void setCurrentPrefetch(int prefetch);

    void prefetchLogLikelihoods();

    void acceptPrefetch(int prefetch);
}
