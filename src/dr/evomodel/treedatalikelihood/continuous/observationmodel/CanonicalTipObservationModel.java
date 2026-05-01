package dr.evomodel.treedatalikelihood.continuous.observationmodel;

import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalTransitionMomentProvider;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;

public interface CanonicalTipObservationModel {

    int getLatentDimension();

    int getObservationDimension();

    TipObservationMode getMode();

    boolean isEmpty();

    /**
     * Fills {@code out} with the upward parent-side canonical message produced by this tip
     * observation model after marginalising or conditioning on the observed tip data through
     * the supplied branch transition.
     *
     * <p>Implementations must write a valid (possibly zero) canonical state into {@code out}
     * and must not modify any other element of the workspace beyond what they own.
     *
     * <p>The caller is responsible for ensuring that {@code gaussianWorkspace} has been
     * allocated for the correct dimension and is not concurrently in use.
     *
     * @param transition       canonical branch transition (child → parent block form)
     * @param momentProvider   moment-form transition provider; only needed for
     *                         {@link TipObservationMode#PARTIAL_EXACT_IDENTITY} — may be
     *                         {@code null} for all other modes
     * @param branchLength     effective branch length; needed when {@code momentProvider}
     *                         is used to fill moment-form matrices on demand
     * @param workspace        tip-parent scratch buffers owned by the traversal
     * @param gaussianWorkspace Gaussian message-ops workspace shared with the traversal
     * @param out              output parent-side canonical state (overwritten)
     */
    void fillParentMessage(CanonicalGaussianTransition transition,
                           CanonicalTransitionMomentProvider momentProvider,
                           double branchLength,
                           TipParentMessageWorkspace workspace,
                           CanonicalGaussianMessageOps.Workspace gaussianWorkspace,
                           CanonicalGaussianState out);

    /**
     * Fills {@code out} with the canonical-form factor for the latent child state induced by
     * this observation model.
     *
     * <p>This is only well-defined for models with finite-precision Gaussian observations
     * (i.e. {@link TipObservationMode#GAUSSIAN_LINK}). Identity-observation models produce
     * degenerate (infinite-precision) factors that cannot be represented as a proper canonical
     * Gaussian state; they must use {@link #fillParentMessage} instead.
     *
     * @throws UnsupportedOperationException if the model does not support a finite canonical
     *                                        child-state factor
     */
    void fillChildCanonicalState(CanonicalGaussianState out,
                                 TipObservationModelWorkspace workspace);

    CanonicalTipObservationModel copy();
}
