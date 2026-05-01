package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianMatrixOps;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.OrthogonalMatrixProvider;
import dr.inference.model.Parameter;
import dr.inference.model.TransposedMatrixParameter;

/**
 * Projection helpers for converting dense canonical selection gradients to the
 * BEAST parameter layout requested by gradient callers.
 */
final class CanonicalSelectionGradientProjector {

    private CanonicalSelectionGradientProjector() { }

    static double[] assembleBlockGradientResult(
            final int dimension,
            final Parameter requestedParameter,
            final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
            final double[] nativeBlockGradient,
            final double[][] gradientR) {
        if (requestedParameter == blockParameter.getParameter()) {
            if (blockParameter.getRotationMatrixParameter() instanceof OrthogonalMatrixProvider) {
                final double[] angleGradient =
                        ((OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter()).pullBackGradient(gradientR);
                final double[] out = new double[blockParameter.getBlockDiagonalNParameters() + angleGradient.length];
                System.arraycopy(nativeBlockGradient, 0, out, 0, blockParameter.getBlockDiagonalNParameters());
                System.arraycopy(angleGradient, 0, out, blockParameter.getBlockDiagonalNParameters(), angleGradient.length);
                return out;
            }
            final double[] out = new double[blockParameter.getBlockDiagonalNParameters() + dimension * dimension];
            System.arraycopy(nativeBlockGradient, 0, out, 0, blockParameter.getBlockDiagonalNParameters());
            final double[] flatRotation = flattenRotationGradient(blockParameter.getRotationMatrixParameter(), gradientR);
            System.arraycopy(flatRotation, 0, out, blockParameter.getBlockDiagonalNParameters(), flatRotation.length);
            return out;
        }

        if (requestedParameter == blockParameter.getRotationMatrixParameter()) {
            return flattenRotationGradient(blockParameter.getRotationMatrixParameter(), gradientR);
        }
        if (blockParameter.getRotationMatrixParameter() instanceof OrthogonalMatrixProvider
                && requestedParameter == ((OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter()).getOrthogonalParameter()) {
            return ((OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter()).pullBackGradient(gradientR);
        }
        if (requestedParameter == blockParameter.getScalarBlockParameter()) {
            return new double[]{nativeBlockGradient[0]};
        }

        final int blockBase = blockParameter.hasLeadingOneByOneBlock() ? 1 : 0;
        final int blockWidth = blockParameter.getNum2x2Blocks();
        for (int family = 0; family < blockParameter.getTwoByTwoParameterFamilyCount(); ++family) {
            if (requestedParameter == blockParameter.getTwoByTwoBlockParameter(family)) {
                final double[] out = new double[blockWidth];
                System.arraycopy(nativeBlockGradient, blockBase + family * blockWidth, out, 0, blockWidth);
                return out;
            }
        }
        throw new IllegalArgumentException("Unsupported block parameter: " + requestedParameter.getId());
    }

    static double[] assembleBlockGradientResultFlat(
            final int dimension,
            final Parameter requestedParameter,
            final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
            final double[] nativeBlockGradient,
            final double[] gradientR) {
        if (requestedParameter == blockParameter.getParameter()) {
            if (blockParameter.getRotationMatrixParameter() instanceof OrthogonalMatrixProvider) {
                final double[] angleGradient =
                        ((OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter())
                                .pullBackGradientFlat(gradientR, dimension);
                final double[] out = new double[blockParameter.getBlockDiagonalNParameters() + angleGradient.length];
                System.arraycopy(nativeBlockGradient, 0, out, 0, blockParameter.getBlockDiagonalNParameters());
                System.arraycopy(angleGradient, 0, out, blockParameter.getBlockDiagonalNParameters(), angleGradient.length);
                return out;
            }
            final double[] out = new double[blockParameter.getBlockDiagonalNParameters() + dimension * dimension];
            System.arraycopy(nativeBlockGradient, 0, out, 0, blockParameter.getBlockDiagonalNParameters());
            final double[] flatRotation = flattenRotationGradientFlat(
                    blockParameter.getRotationMatrixParameter(), gradientR, dimension);
            System.arraycopy(flatRotation, 0, out, blockParameter.getBlockDiagonalNParameters(), flatRotation.length);
            return out;
        }

        if (requestedParameter == blockParameter.getRotationMatrixParameter()) {
            return flattenRotationGradientFlat(blockParameter.getRotationMatrixParameter(), gradientR, dimension);
        }
        if (blockParameter.getRotationMatrixParameter() instanceof OrthogonalMatrixProvider
                && requestedParameter == ((OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter()).getOrthogonalParameter()) {
            return ((OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter())
                    .pullBackGradientFlat(gradientR, dimension);
        }
        if (requestedParameter == blockParameter.getScalarBlockParameter()) {
            return new double[]{nativeBlockGradient[0]};
        }

        final int blockBase = blockParameter.hasLeadingOneByOneBlock() ? 1 : 0;
        final int blockWidth = blockParameter.getNum2x2Blocks();
        for (int family = 0; family < blockParameter.getTwoByTwoParameterFamilyCount(); ++family) {
            if (requestedParameter == blockParameter.getTwoByTwoBlockParameter(family)) {
                final double[] out = new double[blockWidth];
                System.arraycopy(nativeBlockGradient, blockBase + family * blockWidth, out, 0, blockWidth);
                return out;
            }
        }
        throw new IllegalArgumentException("Unsupported block parameter: " + requestedParameter.getId());
    }

    static double[] reorderDenseGradientForRequestedParameter(final int dimension,
                                                              final Parameter requestedParameter,
                                                              final double[] rowMajorGradient) {
        if (requestedParameter instanceof TransposedMatrixParameter) {
            return rowMajorGradient;
        }
        if (!(requestedParameter instanceof dr.inference.model.MatrixParameterInterface)) {
            return rowMajorGradient;
        }
        final double[] reordered = new double[rowMajorGradient.length];
        int index = 0;
        for (int col = 0; col < dimension; ++col) {
            for (int row = 0; row < dimension; ++row) {
                reordered[index++] = rowMajorGradient[row * dimension + col];
            }
        }
        return reordered;
    }

    static double[] pullBackDenseGradientToBlock(
            final int dimension,
            final Parameter requestedParameter,
            final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
            final double[] denseGradient) {
        final int d = dimension;
        final double[][] gradientA = new double[d][d];
        final double[] rData = new double[d * d];
        final double[] rinvData = new double[d * d];
        final double[] rawBlockSource = new double[blockParameter.getCompressedDDimension()];
        final double[] nativeBlockGradient = new double[blockParameter.getBlockDiagonalNParameters()];
        final double[][] gradientD = new double[d][d];
        final double[][] gradientR = new double[d][d];
        final double[][] aMatrix = blockParameter.getParameterAsMatrix();

        fillMatrixFromRowMajorGradient(gradientA, denseGradient, d);
        blockParameter.fillRAndRinv(rData, rinvData);
        computeGradientWrtBlockDiagonalBasis(gradientA, rData, rinvData, gradientD);
        computeGradientWrtRotationMatrix(gradientA, aMatrix, rData, rinvData, gradientR);
        compressActiveBlockGradient(blockParameter, gradientD, rawBlockSource);
        blockParameter.chainGradient(rawBlockSource, nativeBlockGradient);
        return assembleBlockGradientResult(dimension, requestedParameter, blockParameter, nativeBlockGradient, gradientR);
    }

    private static double[] flattenColumnMajor(final double[][] matrix) {
        final int d = matrix.length;
        final double[] out = new double[d * d];
        GaussianMatrixOps.matrixToColumnMajorParameter(matrix, out, d);
        return out;
    }

    private static double[] flattenRowMajor(final double[][] matrix) {
        final int d = matrix.length;
        final double[] out = new double[d * d];
        GaussianMatrixOps.matrixToRowMajor(matrix, out, d);
        return out;
    }

    private static double[] flattenRotationGradient(final Parameter rotationParameter,
                                                    final double[][] gradientR) {
        if (rotationParameter instanceof TransposedMatrixParameter) {
            return flattenRowMajor(gradientR);
        }
        return flattenColumnMajor(gradientR);
    }

    private static double[] flattenRotationGradientFlat(final Parameter rotationParameter,
                                                        final double[] gradientR,
                                                        final int dimension) {
        final double[] out = new double[dimension * dimension];
        if (rotationParameter instanceof TransposedMatrixParameter) {
            System.arraycopy(gradientR, 0, out, 0, out.length);
        } else {
            for (int row = 0; row < dimension; ++row) {
                for (int col = 0; col < dimension; ++col) {
                    out[col * dimension + row] = gradientR[row * dimension + col];
                }
            }
        }
        return out;
    }

    private static void fillMatrixFromRowMajorGradient(final double[][] out,
                                                       final double[] gradient,
                                                       final int d) {
        for (int row = 0; row < d; ++row) {
            for (int col = 0; col < d; ++col) {
                out[row][col] = gradient[row * d + col];
            }
        }
    }

    private static void computeGradientWrtBlockDiagonalBasis(final double[][] gradientA,
                                                             final double[] rData,
                                                             final double[] rinvData,
                                                             final double[][] out) {
        final int d = gradientA.length;
        final double[][] temp = new double[d][d];
        final double[][] rinvTranspose = new double[d][d];

        for (int row = 0; row < d; ++row) {
            for (int col = 0; col < d; ++col) {
                rinvTranspose[row][col] = rinvData[col * d + row];
            }
        }

        multiplyTransposeLeft(rData, d, gradientA, temp);
        multiplyDense(temp, rinvTranspose, out);
    }

    private static void computeGradientWrtRotationMatrix(final double[][] gradientA,
                                                         final double[][] aMatrix,
                                                         final double[] rData,
                                                         final double[] rinvData,
                                                         final double[][] out) {
        final int d = gradientA.length;
        final double[][] aTranspose = new double[d][d];
        final double[][] dMatrix = new double[d][d];
        final double[][] dTranspose = new double[d][d];
        final double[][] temp = new double[d][d];
        final double[][] temp2 = new double[d][d];
        final double[][] rinv = new double[d][d];
        final double[][] rinvTranspose = new double[d][d];
        final double[][] rMatrix = new double[d][d];

        transpose(aMatrix, aTranspose);

        for (int row = 0; row < d; ++row) {
            for (int col = 0; col < d; ++col) {
                rinv[row][col] = rinvData[row * d + col];
                rinvTranspose[row][col] = rinvData[col * d + row];
                rMatrix[row][col] = rData[row * d + col];
            }
        }

        multiplyDense(rinv, aMatrix, temp2);
        multiplyDense(temp2, rMatrix, dMatrix);
        transpose(dMatrix, dTranspose);

        multiplyDense(gradientA, rinvTranspose, temp);
        multiplyDense(temp, dTranspose, temp2);

        multiplyDense(aTranspose, gradientA, temp);
        multiplyDense(temp, rinvTranspose, out);

        subtractIntoOut(temp2, out);
    }

    private static void compressActiveBlockGradient(final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                                    final double[][] gradientD,
                                                    final double[] out) {
        final int d = gradientD.length;
        int blockIndex = 0;
        final int upperBase = d;
        final int lowerBase = d + blockParameter.getNum2x2Blocks();

        for (int i = 0; i < d; ++i) {
            out[i] = gradientD[i][i];
        }

        for (int b = 0; b < blockParameter.getNumBlocks(); ++b) {
            if (blockParameter.getBlockSizes()[b] != 2) {
                continue;
            }
            final int start = blockParameter.getBlockStarts()[b];
            out[upperBase + blockIndex] = gradientD[start][start + 1];
            out[lowerBase + blockIndex] = gradientD[start + 1][start];
            blockIndex++;
        }
    }

    private static void multiplyTransposeLeft(final double[] leftRowMajor,
                                              final int dimension,
                                              final double[][] right,
                                              final double[][] out) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += leftRowMajor[k * dimension + i] * right[k][j];
                }
                out[i][j] = sum;
            }
        }
    }

    private static void multiplyDense(final double[][] left,
                                      final double[][] right,
                                      final double[][] out) {
        final int rows = left.length;
        final int cols = right[0].length;
        final int inner = right.length;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double sum = 0.0;
                for (int k = 0; k < inner; ++k) {
                    sum += left[i][k] * right[k][j];
                }
                out[i][j] = sum;
            }
        }
    }

    private static void transpose(final double[][] in,
                                  final double[][] out) {
        for (int i = 0; i < in.length; ++i) {
            for (int j = 0; j < in[i].length; ++j) {
                out[j][i] = in[i][j];
            }
        }
    }

    private static void subtractIntoOut(final double[][] left,
                                        final double[][] out) {
        for (int i = 0; i < left.length; ++i) {
            for (int j = 0; j < left[i].length; ++j) {
                out[i][j] = left[i][j] - out[i][j];
            }
        }
    }
}
