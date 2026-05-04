package dr.evomodel.treedatalikelihood.continuous.canonical.gradient;

import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.inference.model.Parameter;

public final class DenseCanonicalSelectionPullback implements CanonicalSelectionPullback {

    private final int dimension;
    private final Parameter requestedParameter;

    public DenseCanonicalSelectionPullback(final int dimension,
                                           final Parameter requestedParameter) {
        this.dimension = dimension;
        this.requestedParameter = requestedParameter;
    }

    @Override
    public double[] gradientForBranch(final double branchLength,
                                      final double[] optimum,
                                      final BranchSufficientStatistics statistics,
                                      final DenseSelectionGradientProvider denseGradientProvider) {
        return projectDenseGradient(denseGradientProvider.gradient(branchLength, optimum, statistics));
    }

    @Override
    public double[] projectDenseGradient(final double[] denseGradient) {
        if (requestedParameter == null) {
            return denseGradient.clone();
        }
        return CanonicalSelectionGradientProjector.reorderDenseGradientForRequestedParameter(
                dimension, requestedParameter, denseGradient);
    }
}
