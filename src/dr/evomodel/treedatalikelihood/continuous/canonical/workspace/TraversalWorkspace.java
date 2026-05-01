package dr.evomodel.treedatalikelihood.continuous.canonical.workspace;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.TipObservationModelWorkspace;

public final class TraversalWorkspace {
    public final double[][] traitCovariance;
    public final CanonicalGaussianTransition transition;
    public final CanonicalGaussianState state;
    public final CanonicalGaussianState siblingProduct;
    public final CanonicalGaussianState downwardParentState;
    public final CanonicalGaussianMessageOps.Workspace gaussianWorkspace;
    public final TipObservationModelWorkspace observationWorkspace;

    public TraversalWorkspace(final int dim) {
        this.traitCovariance = new double[dim][dim];
        this.transition = new CanonicalGaussianTransition(dim);
        this.state = new CanonicalGaussianState(dim);
        this.siblingProduct = new CanonicalGaussianState(dim);
        this.downwardParentState = new CanonicalGaussianState(dim);
        this.gaussianWorkspace = new CanonicalGaussianMessageOps.Workspace(dim);
        this.observationWorkspace = new TipObservationModelWorkspace();
    }
}
