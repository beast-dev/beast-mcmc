package dr.evomodel.treedatalikelihood.continuous.integration;

import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalOUTransitionProvider;

final class CanonicalOUProviderSupport {

    private CanonicalOUProviderSupport() {
        // no instances
    }

    static CanonicalOUTransitionProvider requireOUProvider(
            final CanonicalBranchTransitionProvider transitionProvider) {
        if (!(transitionProvider instanceof CanonicalOUTransitionProvider)) {
            throw new UnsupportedOperationException(
                    "Not yet implemented: canonical OU gradients currently support only "
                            + "CanonicalOUTransitionProvider implementations.");
        }
        return (CanonicalOUTransitionProvider) transitionProvider;
    }
}
