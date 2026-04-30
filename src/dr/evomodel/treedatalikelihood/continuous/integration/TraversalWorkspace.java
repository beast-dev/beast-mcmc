package dr.evomodel.treedatalikelihood.continuous.integration;

import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalGaussianMessageOps;

final class TraversalWorkspace {
    final double[][] traitCovariance;
    final CanonicalGaussianTransition transition;
    final CanonicalGaussianState state;
    final CanonicalGaussianState siblingProduct;
    final CanonicalGaussianState downwardParentState;
    final CanonicalGaussianMessageOps.Workspace gaussianWorkspace;

    TraversalWorkspace(final int dim) {
        this.traitCovariance = new double[dim][dim];
        this.transition = new CanonicalGaussianTransition(dim);
        this.state = new CanonicalGaussianState(dim);
        this.siblingProduct = new CanonicalGaussianState(dim);
        this.downwardParentState = new CanonicalGaussianState(dim);
        this.gaussianWorkspace = new CanonicalGaussianMessageOps.Workspace(dim);
    }
}
