package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.Parameter;

final class NativeSpecializedCanonicalSelectionPullback implements CanonicalSelectionPullback {

    private final OUCanonicalNativeSelectionGradient nativeSelectionGradient;
    private final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter;
    private final Parameter requestedParameter;

    NativeSpecializedCanonicalSelectionPullback(
            final OUCanonicalNativeSelectionGradient nativeSelectionGradient,
            final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
            final Parameter requestedParameter) {
        this.nativeSelectionGradient = nativeSelectionGradient;
        this.blockParameter = blockParameter;
        this.requestedParameter = requestedParameter;
    }

    @Override
    public double[] gradientForBranch(final double branchLength,
                                      final double[] optimum,
                                      final BranchSufficientStatistics statistics,
                                      final DenseSelectionGradientProvider denseGradientProvider) {
        return nativeSelectionGradient.getNativeSpecializedGradientWrtSelection(
                branchLength,
                optimum,
                statistics,
                blockParameter,
                requestedParameter);
    }

    @Override
    public double[] projectDenseGradient(final double[] denseGradient) {
        throw new UnsupportedOperationException(
                "Native specialized selection pullback needs branch-local canonical adjoints.");
    }
}
