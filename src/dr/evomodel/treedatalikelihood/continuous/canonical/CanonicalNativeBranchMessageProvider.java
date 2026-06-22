package dr.evomodel.treedatalikelihood.continuous.canonical;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.BranchGradientWorkspace;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.CanonicalTipObservationModel;

/**
 * Optional provider capability for branch models that can propagate canonical messages
 * without first materializing a full two-block canonical transition object.
 */
public interface CanonicalNativeBranchMessageProvider {

    boolean fillBackwardMessageForState(int childNodeIndex,
                                        CanonicalGaussianState childState,
                                        BranchGradientWorkspace workspace,
                                        CanonicalGaussianState out);

    boolean fillBackwardMessageForTip(int childNodeIndex,
                                      CanonicalTipObservationModel observationModel,
                                      BranchGradientWorkspace workspace,
                                      CanonicalGaussianState out);

    boolean fillForwardMessage(int childNodeIndex,
                               CanonicalGaussianState parentState,
                               BranchGradientWorkspace workspace,
                               CanonicalGaussianState out);

    default boolean bypassCanonicalTransitionPreload() {
        return true;
    }
}
