package dr.evomodel.treedatalikelihood.continuous.canonical;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.TipObservationModelWorkspace;

final class TraversalWorkspace {
    final double[][] traitCovariance;
    final CanonicalGaussianTransition transition;
    final CanonicalGaussianState state;
    final CanonicalGaussianState siblingProduct;
    final CanonicalGaussianState downwardParentState;
    final CanonicalGaussianMessageOps.Workspace gaussianWorkspace;
    final TipObservationModelWorkspace observationWorkspace;

    TraversalWorkspace(final int dim) {
        this.traitCovariance = new double[dim][dim];
        this.transition = new CanonicalGaussianTransition(dim);
        this.state = new CanonicalGaussianState(dim);
        this.siblingProduct = new CanonicalGaussianState(dim);
        this.downwardParentState = new CanonicalGaussianState(dim);
        this.gaussianWorkspace = new CanonicalGaussianMessageOps.Workspace(dim);
        this.observationWorkspace = new TipObservationModelWorkspace();
    }
}
