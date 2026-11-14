package dr.evomodel.treedatalikelihood.continuousNew.observationmodel;

import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianState;

/**
 * Narrow branch-message view needed by tip observation models.
 */
public interface CanonicalTipBranchMessageProvider {

    boolean fillBackwardMessageForObservedTip(int childNodeIndex,
                                              double[] observedValue,
                                              CanonicalGaussianState out);
}
