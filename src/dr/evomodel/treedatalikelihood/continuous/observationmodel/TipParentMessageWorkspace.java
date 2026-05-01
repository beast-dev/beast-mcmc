package dr.evomodel.treedatalikelihood.continuous.observationmodel;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;

/**
 * Shared scratch buffers for computing the upward parent message from a tip observation model.
 *
 * <p>Owned by the traversal workspace and passed to
 * {@link CanonicalTipObservationModel#fillParentMessage}. Each implementation uses only the
 * fields it needs; unused fields are ignored.
 *
 * <ul>
 *   <li>{@link #childStateScratch} — intermediate canonical state for Gaussian-link models
 *       (filled by {@code fillChildCanonicalState} then pushed backward).</li>
 *   <li>{@link #observationModelWorkspace} — covariance inversion scratch for Gaussian-link models.</li>
 *   <li>{@link #partialIdentityProjection} — Schur-complement projection for partially observed
 *       exact identity tips.</li>
 * </ul>
 *
 * <p>Note: the Gaussian message-passing workspace ({@code CanonicalGaussianMessageOps.Workspace})
 * is shared with the rest of the traversal and is passed separately to {@code fillParentMessage}
 * rather than being stored here.
 */
public final class TipParentMessageWorkspace {

    public final CanonicalGaussianState childStateScratch;
    public final TipObservationModelWorkspace observationModelWorkspace;
    public final PartialIdentityTipProjection partialIdentityProjection;

    public TipParentMessageWorkspace(final int dim) {
        if (dim < 1) throw new IllegalArgumentException("dim must be >= 1");
        this.childStateScratch         = new CanonicalGaussianState(dim);
        this.observationModelWorkspace  = new TipObservationModelWorkspace();
        this.partialIdentityProjection  = new PartialIdentityTipProjection(dim);
    }
}
