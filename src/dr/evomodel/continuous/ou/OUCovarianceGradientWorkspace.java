package dr.evomodel.continuous.ou;

final class OUCovarianceGradientWorkspace {

    final int dimension;
    final int squareStride;
    final int blockDimension;
    final int blockStride;
    final int lyapunovAugmentedRows;
    final int lyapunovAugmentedColumns;
    final double[] squareMatrices;
    final double[] blockMatrices;
    final double[] lyapunovAugmented;

    OUCovarianceGradientWorkspace(final int dimension,
                                  final int squareMatrixCount,
                                  final int blockMatrixCount,
                                  final boolean allocateLyapunovAugmented) {
        this.dimension = dimension;
        this.squareStride = dimension * dimension;
        this.blockDimension = 2 * dimension;
        this.blockStride = blockDimension * blockDimension;
        this.squareMatrices = OUCovarianceGradientMath.allocateMatrixStack(
                squareMatrixCount, squareStride);
        this.blockMatrices = OUCovarianceGradientMath.allocateMatrixStack(
                blockMatrixCount, blockStride);

        final int lyapunovDimension = dimension * dimension;
        this.lyapunovAugmentedRows = lyapunovDimension;
        this.lyapunovAugmentedColumns = lyapunovDimension + 1;
        this.lyapunovAugmented = allocateLyapunovAugmented
                ? OUCovarianceGradientMath.allocateMatrix(lyapunovAugmentedRows, lyapunovAugmentedColumns)
                : null;
    }

    int squareOffset(final int index) {
        return index * squareStride;
    }

    int blockOffset(final int index) {
        return index * blockStride;
    }
}
