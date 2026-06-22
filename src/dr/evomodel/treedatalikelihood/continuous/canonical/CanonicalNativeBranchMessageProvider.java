package dr.evomodel.treedatalikelihood.continuous.canonical;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.BranchGradientWorkspace;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.CanonicalTipBranchMessageProvider;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.TipObservationMode;

/**
 * Optional provider capability for branch models that can propagate canonical messages
 * without first materializing a full two-block canonical transition object.
 */
public interface CanonicalNativeBranchMessageProvider extends CanonicalTipBranchMessageProvider {

    boolean fillBackwardMessageForState(int childNodeIndex,
                                        CanonicalGaussianState childState,
                                        BranchGradientWorkspace workspace,
                                        CanonicalGaussianState out);

    boolean fillForwardMessage(int childNodeIndex,
                               CanonicalGaussianState parentState,
                               BranchGradientWorkspace workspace,
                               CanonicalGaussianState out);

    default boolean hasNativeBranchMessages() {
        return true;
    }

    default boolean supportsPostorderInternalMessage() {
        return hasNativeBranchMessages();
    }

    default boolean supportsPostorderTipMessage(final TipObservationMode mode) {
        return false;
    }

    default boolean supportsPreorderForwardMessage() {
        return false;
    }

    default boolean supportsFixedRootForwardMessage() {
        return false;
    }

    default void clearBranchMessageCache() {
    }
}
