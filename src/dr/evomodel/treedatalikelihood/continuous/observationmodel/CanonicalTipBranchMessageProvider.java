package dr.evomodel.treedatalikelihood.continuous.observationmodel;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;

/**
 * Narrow branch-message view needed by tip observation models.
 */
public interface CanonicalTipBranchMessageProvider {

    boolean fillBackwardMessageForObservedTip(int childNodeIndex,
                                              double[] observedValue,
                                              CanonicalGaussianState out);
}
