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
            final double[] angleGradientScratch,
            final double[] gradientAccumulator) {
        final int nativeBlockDim = blockParameter.getBlockDiagonalNParameters();
        blockParameter.chainGradient(compressedBlockGradient, nativeBlockGradientScratch);
        final int angleDim = orthogonalRotation.getOrthogonalParameter().getDimension();
        orthogonalRotation.fillPullBackGradientFlat(
                rotationGradient,
                blockParameter.getRowDimension(),
                angleGradientScratch);
        for (int i = 0; i < nativeBlockDim; ++i) {
            gradientAccumulator[i] += nativeBlockGradientScratch[i];
        }
        for (int i = 0; i < angleDim; ++i) {
            gradientAccumulator[nativeBlockDim + i] += angleGradientScratch[i];
        }
    }
}
