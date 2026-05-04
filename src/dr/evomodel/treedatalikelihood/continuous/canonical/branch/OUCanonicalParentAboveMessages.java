package dr.evomodel.treedatalikelihood.continuous.canonical.branch;

import dr.evomodel.treedatalikelihood.continuous.CanonicalDebugOptions;
import dr.evomodel.treedatalikelihood.continuous.CanonicalGradientFallbackPolicy;
import dr.evomodel.treedatalikelihood.continuous.OUGaussianBranchTransitionProvider;

import dr.evomodel.treedatalikelihood.preorder.NormalSufficientStatistics;
import org.ejml.data.DenseMatrix64F;

final class OUCanonicalParentAboveMessages {

    private final OUCanonicalTransitionState transitionState;

    OUCanonicalParentAboveMessages(final OUCanonicalTransitionState transitionState) {
        this.transitionState = transitionState;
    }

    DenseMatrix64F recoverOrRequire(final NormalSufficientStatistics aboveChild,
                                    final NormalSufficientStatistics aboveParent) {
        return transitionState.recoverOrUseParentAbovePrecision(aboveChild, aboveParent);
    }

    DenseMatrix64F require(final NormalSufficientStatistics aboveParent) {
        return transitionState.requireParentAbovePrecision(aboveParent);
    }
}
