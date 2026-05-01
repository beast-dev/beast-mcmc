package dr.evomodel.treedatalikelihood.continuous.canonical;

import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalOUTransitionProvider;

public final class CanonicalOUProviderSupport {

    private CanonicalOUProviderSupport() {
        // no instances
    }

    public static CanonicalOUTransitionProvider requireOUProvider(
            final CanonicalBranchTransitionProvider transitionProvider) {
        if (!(transitionProvider instanceof CanonicalOUTransitionProvider)) {
            throw new UnsupportedOperationException(
                    "Not yet implemented: canonical OU gradients currently support only "
                            + "CanonicalOUTransitionProvider implementations.");
        }
        return (CanonicalOUTransitionProvider) transitionProvider;
    }
}
