package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.preorder.NormalSufficientStatistics;
import org.ejml.data.DenseMatrix64F;

final class OUCanonicalParentAboveMessages {

    private final OUCanonicalTransitionState transitionState;

    OUCanonicalParentAboveMessages(final OUCanonicalTransitionState transitionState) {
        this.transitionState = transitionState;
    }

    DenseMatrix64F require(final NormalSufficientStatistics aboveParent) {
        return transitionState.requireParentAbovePrecision(aboveParent);
    }
}
