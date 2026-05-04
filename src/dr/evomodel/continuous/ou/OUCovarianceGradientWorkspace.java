package dr.evomodel.continuous.ou;

final class OUCovarianceGradientWorkspace {

    final int dimension;
    final double[][][] squareMatrices;
    final double[][][] blockMatrices;
    final double[][] lyapunovAugmented;

    OUCovarianceGradientWorkspace(final int dimension,
                                  final int squareMatrixCount,
                                  final int blockMatrixCount,
                                  final boolean allocateLyapunovAugmented) {
        this.dimension = dimension;
        this.squareMatrices = OUCovarianceGradientMath.allocateMatrixStack(
                squareMatrixCount, dimension, dimension);
        this.blockMatrices = OUCovarianceGradientMath.allocateMatrixStack(
                blockMatrixCount, 2 * dimension, 2 * dimension);
        this.lyapunovAugmented = allocateLyapunovAugmented
                ? OUCovarianceGradientMath.allocateMatrix(dimension * dimension, dimension * dimension + 1)
                : null;
    }
}
