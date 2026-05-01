package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.OrthogonalMatrixProvider;

final class OrthogonalBlockGradientPullback {

    private OrthogonalBlockGradientPullback() { }

    static void accumulateNativeParameterGradient(
            final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
            final OrthogonalMatrixProvider orthogonalRotation,
            final double[] compressedBlockGradient,
            final double[] rotationGradient,
            final double[] nativeBlockGradientScratch,
            final double[] gradientAccumulator) {
        final int nativeBlockDim = blockParameter.getBlockDiagonalNParameters();
        blockParameter.chainGradient(compressedBlockGradient, nativeBlockGradientScratch);
        final double[] angleGradient =
                orthogonalRotation.pullBackGradientFlat(rotationGradient, blockParameter.getRowDimension());
        for (int i = 0; i < nativeBlockDim; ++i) {
            gradientAccumulator[i] += nativeBlockGradientScratch[i];
        }
        for (int i = 0; i < angleGradient.length; ++i) {
            gradientAccumulator[nativeBlockDim + i] += angleGradient[i];
        }
    }
}
