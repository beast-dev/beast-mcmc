package dr.evomodel.treedatalikelihood.continuousNew.observationmodel;

import dr.evomodel.continuous.canonical.workspace.BranchGradientWorkspace;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalTransitionMomentProvider;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianState;

public final class LinearExactTipContribution {

    private final int dimension;

    public LinearExactTipContribution(final int dimension) {
        if (dimension < 1) {
            throw new IllegalArgumentException("dimension must be positive");
        }
        this.dimension = dimension;
    }

    public void fillLocalAdjointsForLinearTip(
            final CanonicalGaussianState aboveState,
            final LinearExactTipObservationModel tipObservation,
            final int childNodeIndex,
            final CanonicalTransitionMomentProvider momentProvider,
            final double branchLength,
            final BranchGradientWorkspace workspace) {
        workspace.linearExactProjection.projectObservedLinearChildToParent(
                tipObservation,
                childNodeIndex,
                momentProvider,
                branchLength,
                workspace.state);
        CanonicalGaussianMessageOps.combineStates(
                aboveState,
                workspace.state,
                workspace.parentPosterior);
        workspace.adjoint.fillMomentsFromCanonical(
                workspace.parentPosterior,
                workspace.mean,
                workspace.covariance,
                dimension);
        fillSecondMomentPrecisionSeed(workspace.mean, workspace.covariance, workspace.varianceFlat);
        workspace.linearExactProjection.fillLocalAdjointsForCurrentProjection(
                tipObservation,
                workspace.mean,
                workspace.varianceFlat,
                -1.0,
                workspace.adjoints);
    }

    public void fillLocalAdjointsForFixedParentLinearTip(
            final LinearExactTipObservationModel tipObservation,
            final int childNodeIndex,
            final CanonicalTransitionMomentProvider momentProvider,
            final double branchLength,
            final double[] fixedRootValue,
            final BranchGradientWorkspace workspace) {
        fillFixedParentPrecisionSeed(fixedRootValue, workspace.varianceFlat);
        workspace.linearExactProjection.fillLocalAdjoints(
                tipObservation,
                childNodeIndex,
                momentProvider,
                branchLength,
                fixedRootValue,
                workspace.varianceFlat,
                -1.0,
                workspace.adjoints);
    }

    private void fillSecondMomentPrecisionSeed(final double[] mean,
                                               final double[] covariance,
                                               final double[] out) {
        for (int i = 0; i < dimension; ++i) {
            final double meanI = mean[i];
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                out[rowOffset + j] = -0.5 * (covariance[rowOffset + j] + meanI * mean[j]);
            }
        }
    }

    private void fillFixedParentPrecisionSeed(final double[] fixedRootValue,
                                              final double[] out) {
        for (int i = 0; i < dimension; ++i) {
            final double xi = fixedRootValue[i];
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                out[rowOffset + j] = -0.5 * xi * fixedRootValue[j];
            }
        }
    }
}
