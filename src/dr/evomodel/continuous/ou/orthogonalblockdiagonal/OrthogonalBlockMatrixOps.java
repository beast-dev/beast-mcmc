package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

final class OrthogonalBlockMatrixOps {

    private OrthogonalBlockMatrixOps() {
        // no instances
    }

    static void multiplyRightBlockDiagonal(final double[] matrixData,
                                           final double[] compressedBlockData,
                                           final int dimension,
                                           final int[] blockStarts,
                                           final int[] blockSizes,
                                           final double[] outData) {
        int block2x2Index = 0;
        final int upperBase = dimension;
        final int lowerBase = dimension + countTwoByTwoBlocks(blockSizes);
        for (int row = 0; row < dimension; ++row) {
            final int rowOffset = row * dimension;
            for (int b = 0; b < blockStarts.length; ++b) {
                final int start = blockStarts[b];
                final int size = blockSizes[b];
                if (size == 1) {
                    outData[rowOffset + start] = matrixData[rowOffset + start] * compressedBlockData[start];
                } else {
                    final double e00 = compressedBlockData[start];
                    final double e11 = compressedBlockData[start + 1];
                    final double e01 = compressedBlockData[upperBase + block2x2Index];
                    final double e10 = compressedBlockData[lowerBase + block2x2Index];
                    final double x0 = matrixData[rowOffset + start];
                    final double x1 = matrixData[rowOffset + start + 1];
                    outData[rowOffset + start] = x0 * e00 + x1 * e10;
                    outData[rowOffset + start + 1] = x0 * e01 + x1 * e11;
                    block2x2Index++;
                }
            }
            block2x2Index = 0;
        }
    }

    static void multiplyBlockDiagonalLeft(final double[] compressedBlockData,
                                          final double[] matrixData,
                                          final double[] outData,
                                          final int dimension,
                                          final int[] blockStarts,
                                          final int[] blockSizes) {
        int block2x2Index = 0;
        final int upperBase = dimension;
        final int lowerBase = dimension + countTwoByTwoBlocks(blockSizes);
        for (int b = 0; b < blockStarts.length; ++b) {
            final int start = blockStarts[b];
            final int size = blockSizes[b];
            if (size == 1) {
                final int rowOffset = start * dimension;
                for (int col = 0; col < dimension; ++col) {
                    outData[rowOffset + col] = compressedBlockData[start] * matrixData[rowOffset + col];
                }
            } else {
                final int row0 = start * dimension;
                final int row1 = (start + 1) * dimension;
                final double e00 = compressedBlockData[start];
                final double e11 = compressedBlockData[start + 1];
                final double e01 = compressedBlockData[upperBase + block2x2Index];
                final double e10 = compressedBlockData[lowerBase + block2x2Index];
                for (int col = 0; col < dimension; ++col) {
                    final double x0 = matrixData[row0 + col];
                    final double x1 = matrixData[row1 + col];
                    outData[row0 + col] = e00 * x0 + e01 * x1;
                    outData[row1 + col] = e10 * x0 + e11 * x1;
                }
                block2x2Index++;
            }
        }
    }

    static void multiplyRightBlockDiagonalTranspose(final double[] matrixData,
                                                    final double[] compressedBlockData,
                                                    final double[] outData,
                                                    final int dimension,
                                                    final int[] blockStarts,
                                                    final int[] blockSizes) {
        int block2x2Index = 0;
        final int upperBase = dimension;
        final int lowerBase = dimension + countTwoByTwoBlocks(blockSizes);
        for (int row = 0; row < dimension; ++row) {
            final int rowOffset = row * dimension;
            for (int b = 0; b < blockStarts.length; ++b) {
                final int start = blockStarts[b];
                final int size = blockSizes[b];
                if (size == 1) {
                    outData[rowOffset + start] = matrixData[rowOffset + start] * compressedBlockData[start];
                } else {
                    final double e00 = compressedBlockData[start];
                    final double e11 = compressedBlockData[start + 1];
                    final double e01 = compressedBlockData[upperBase + block2x2Index];
                    final double e10 = compressedBlockData[lowerBase + block2x2Index];
                    final double x0 = matrixData[rowOffset + start];
                    final double x1 = matrixData[rowOffset + start + 1];
                    outData[rowOffset + start] = x0 * e00 + x1 * e01;
                    outData[rowOffset + start + 1] = x0 * e10 + x1 * e11;
                    block2x2Index++;
                }
            }
            block2x2Index = 0;
        }
    }

    static void multiplyBlockDiagonalLeftTranspose(final double[] compressedBlockData,
                                                   final double[] matrixData,
                                                   final double[] outData,
                                                   final int dimension,
                                                   final int[] blockStarts,
                                                   final int[] blockSizes) {
        int block2x2Index = 0;
        final int upperBase = dimension;
        final int lowerBase = dimension + countTwoByTwoBlocks(blockSizes);
        for (int b = 0; b < blockStarts.length; ++b) {
            final int start = blockStarts[b];
            final int size = blockSizes[b];
            if (size == 1) {
                final int rowOffset = start * dimension;
                for (int col = 0; col < dimension; ++col) {
                    outData[rowOffset + col] = compressedBlockData[start] * matrixData[rowOffset + col];
                }
            } else {
                final int row0 = start * dimension;
                final int row1 = (start + 1) * dimension;
                final double e00 = compressedBlockData[start];
                final double e11 = compressedBlockData[start + 1];
                final double e01 = compressedBlockData[upperBase + block2x2Index];
                final double e10 = compressedBlockData[lowerBase + block2x2Index];
                for (int col = 0; col < dimension; ++col) {
                    final double x0 = matrixData[row0 + col];
                    final double x1 = matrixData[row1 + col];
                    outData[row0 + col] = e00 * x0 + e10 * x1;
                    outData[row1 + col] = e01 * x0 + e11 * x1;
                }
                block2x2Index++;
            }
        }
    }

    static void applyBlockDiagonalTranspose(final double[] compressedBlockData,
                                            final double[] in,
                                            final double[] out,
                                            final int dimension,
                                            final int[] blockStarts,
                                            final int[] blockSizes) {
        int block2x2Index = 0;
        final int upperBase = dimension;
        final int lowerBase = dimension + countTwoByTwoBlocks(blockSizes);
        for (int b = 0; b < blockStarts.length; ++b) {
            final int start = blockStarts[b];
            final int size = blockSizes[b];
            if (size == 1) {
                out[start] = compressedBlockData[start] * in[start];
            } else {
                final double e00 = compressedBlockData[start];
                final double e11 = compressedBlockData[start + 1];
                final double e01 = compressedBlockData[upperBase + block2x2Index];
                final double e10 = compressedBlockData[lowerBase + block2x2Index];
                final double x0 = in[start];
                final double x1 = in[start + 1];
                out[start] = e00 * x0 + e10 * x1;
                out[start + 1] = e01 * x0 + e11 * x1;
                block2x2Index++;
            }
        }
    }

    private static int countTwoByTwoBlocks(final int[] blockSizes) {
        int count = 0;
        for (int size : blockSizes) {
            if (size == 2) {
                count++;
            }
        }
        return count;
    }
}
