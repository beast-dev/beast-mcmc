package dr.evomodel.treedatalikelihood.continuous.canonical.gradient;

import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.Parameter;

public final class NativeSpecializedCanonicalSelectionPullback implements CanonicalSelectionPullback {

    public interface NativeSelectionGradientProvider {
        double[] getNativeSpecializedGradientWrtSelection(double branchLength,
                                                          double[] optimum,
                                                          BranchSufficientStatistics statistics,
                                                          AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                                          Parameter requestedParameter);
    }

    private final NativeSelectionGradientProvider nativeSelectionGradient;
    private final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter;
    private final Parameter requestedParameter;

    public NativeSpecializedCanonicalSelectionPullback(
            final NativeSelectionGradientProvider nativeSelectionGradient,
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
