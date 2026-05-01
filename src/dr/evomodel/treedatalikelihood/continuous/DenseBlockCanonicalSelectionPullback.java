package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.Parameter;

final class DenseBlockCanonicalSelectionPullback implements CanonicalSelectionPullback {

    private final int dimension;
    private final Parameter requestedParameter;
    private final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter;

    DenseBlockCanonicalSelectionPullback(
            final int dimension,
            final Parameter requestedParameter,
            final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter) {
        this.dimension = dimension;
        this.requestedParameter = requestedParameter;
        this.blockParameter = blockParameter;
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
        return CanonicalSelectionGradientProjector.pullBackDenseGradientToBlock(
                dimension, requestedParameter, blockParameter, denseGradient);
    }
}
