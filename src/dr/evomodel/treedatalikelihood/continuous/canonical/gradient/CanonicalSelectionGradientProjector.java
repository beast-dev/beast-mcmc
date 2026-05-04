package dr.evomodel.treedatalikelihood.continuous.canonical.gradient;

import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.OrthogonalMatrixProvider;
import dr.inference.model.Parameter;
import dr.inference.model.TransposedMatrixParameter;

import java.util.concurrent.atomic.AtomicLong;

/**
 * Projection helpers for converting dense canonical selection gradients to the
 * BEAST parameter layout requested by gradient callers.
 */
public final class CanonicalSelectionGradientProjector {

    private static final AtomicLong LEGACY_MATRIX_PULLBACK_COUNT = new AtomicLong();

    public static final class Workspace {
        private final int dimension;
        private final double[] gradientA;
        private final double[] rData;
        private final double[] rinvData;
        private final double[] rawBlockSource;
        private final double[] nativeBlockGradient;
        private final double[] gradientD;
        private final double[] gradientR;
        private final double[] aMatrix;
        private final double[] temp;
        private final double[] temp2;
        private final double[] temp3;

        public Workspace(final int dimension,
                         final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter) {
            if (dimension < 1) {
                throw new IllegalArgumentException("dimension must be positive");
            }
            this.dimension = dimension;
            final int d2 = dimension * dimension;
            this.gradientA = new double[d2];
            this.rData = new double[d2];
            this.rinvData = new double[d2];
            this.rawBlockSource = new double[blockParameter.getCompressedDDimension()];
            this.nativeBlockGradient = new double[blockParameter.getBlockDiagonalNParameters()];
            this.gradientD = new double[d2];
            this.gradientR = new double[d2];
            this.aMatrix = new double[d2];
            this.temp = new double[d2];
            this.temp2 = new double[d2];
            this.temp3 = new double[d2];
        }

        private void checkCompatible(final int dimension,
                                     final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter) {
            if (this.dimension != dimension) {
                throw new IllegalArgumentException(
                        "workspace dimension must be " + dimension + " but is " + this.dimension);
            }
            if (rawBlockSource.length != blockParameter.getCompressedDDimension()
                    || nativeBlockGradient.length != blockParameter.getBlockDiagonalNParameters()) {
                throw new IllegalArgumentException("workspace is not compatible with the block parameter");
            }
        }
    }

    private CanonicalSelectionGradientProjector() { }

    public static void resetInstrumentation() {
        LEGACY_MATRIX_PULLBACK_COUNT.set(0L);
    }

    public static long getLegacyMatrixPullbackCount() {
        return LEGACY_MATRIX_PULLBACK_COUNT.get();
    }

    public static double[] assembleBlockGradientResultFlat(
            final int dimension,
            final Parameter requestedParameter,
            final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
            final double[] nativeBlockGradient,
            final double[] gradientR) {
        if (requestedParameter == blockParameter.getParameter()) {
            if (blockParameter.getRotationMatrixParameter() instanceof OrthogonalMatrixProvider) {
                final OrthogonalMatrixProvider orthogonalProvider =
                        (OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter();
                final int nativeBlockDim = blockParameter.getBlockDiagonalNParameters();
                final int angleDim = orthogonalProvider.getOrthogonalParameter().getDimension();
                final double[] out = new double[nativeBlockDim + angleDim];
                System.arraycopy(nativeBlockGradient, 0, out, 0, blockParameter.getBlockDiagonalNParameters());
                orthogonalProvider.fillPullBackGradientFlat(gradientR, dimension, out, nativeBlockDim);
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
            final OrthogonalMatrixProvider orthogonalProvider =
                    (OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter();
            final double[] out = new double[orthogonalProvider.getOrthogonalParameter().getDimension()];
            orthogonalProvider.fillPullBackGradientFlat(gradientR, dimension, out);
            return out;
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

    public static double[] reorderDenseGradientForRequestedParameter(final int dimension,
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

    public static double[] pullBackDenseGradientToBlock(
            final int dimension,
            final Parameter requestedParameter,
            final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
            final double[] denseGradient) {
        LEGACY_MATRIX_PULLBACK_COUNT.incrementAndGet();
        return pullBackDenseGradientToBlock(
                dimension,
                requestedParameter,
                blockParameter,
                denseGradient,
                new Workspace(dimension, blockParameter));
    }

    public static double[] pullBackDenseGradientToBlock(
            final int dimension,
            final Parameter requestedParameter,
            final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
            final double[] denseGradient,
            final Workspace workspace) {
        fillBlockGradientProjectionFlat(dimension, blockParameter, denseGradient, workspace);
        return assembleBlockGradientResultFlat(
                dimension,
                requestedParameter,
                blockParameter,
                workspace.nativeBlockGradient,
                workspace.gradientR);
    }

    public static void fillBlockGradientProjectionFlat(
            final int dimension,
            final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
            final double[] denseGradient,
            final Workspace workspace) {
        workspace.checkCompatible(dimension, blockParameter);
        final int d = dimension;
        System.arraycopy(denseGradient, 0, workspace.gradientA, 0, d * d);
        blockParameter.fillRAndRinv(workspace.rData, workspace.rinvData);
        fillBlockParameterMatrixFlat(blockParameter, workspace.aMatrix, d);
        computeGradientWrtBlockDiagonalBasisFlat(
                workspace.gradientA,
                workspace.rData,
                workspace.rinvData,
                workspace.gradientD,
                workspace.temp,
                d);
        computeGradientWrtRotationMatrixFlat(
                workspace.gradientA,
                workspace.aMatrix,
                workspace.rData,
                workspace.rinvData,
                workspace.gradientR,
                workspace.temp,
                workspace.temp2,
                workspace.temp3,
                d);
        compressActiveBlockGradientFlat(blockParameter, workspace.gradientD, workspace.rawBlockSource, d);
        blockParameter.chainGradient(workspace.rawBlockSource, workspace.nativeBlockGradient);
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

    private static void fillBlockParameterMatrixFlat(final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                                     final double[] out,
                                                     final int dimension) {
        for (int row = 0; row < dimension; ++row) {
            final int rowOffset = row * dimension;
            for (int col = 0; col < dimension; ++col) {
                out[rowOffset + col] = blockParameter.getParameterValue(row, col);
            }
        }
    }

    private static void computeGradientWrtBlockDiagonalBasisFlat(final double[] gradientA,
                                                                 final double[] rData,
                                                                 final double[] rinvData,
                                                                 final double[] out,
                                                                 final double[] temp,
                                                                 final int d) {
        multiplyTransposeLeftFlat(rData, gradientA, temp, d);
        multiplyRightTransposeFlat(temp, rinvData, out, d);
    }

    private static void computeGradientWrtRotationMatrixFlat(final double[] gradientA,
                                                            final double[] aMatrix,
                                                            final double[] rData,
                                                            final double[] rinvData,
                                                            final double[] out,
                                                            final double[] temp,
                                                            final double[] temp2,
                                                            final double[] temp3,
                                                            final int d) {
        multiplyFlat(rinvData, aMatrix, temp, d);
        multiplyFlat(temp, rData, temp2, d);

        multiplyRightTransposeFlat(gradientA, rinvData, temp, d);
        multiplyRightTransposeFlat(temp, temp2, temp3, d);

        multiplyTransposeLeftFlat(aMatrix, gradientA, temp, d);
        multiplyRightTransposeFlat(temp, rinvData, out, d);

        for (int i = 0; i < d * d; ++i) {
            out[i] = temp3[i] - out[i];
        }
    }

    private static void compressActiveBlockGradientFlat(final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                                        final double[] gradientD,
                                                        final double[] out,
                                                        final int d) {
        int blockIndex = 0;
        final int upperBase = d;
        final int lowerBase = d + blockParameter.getNum2x2Blocks();

        for (int i = 0; i < d; ++i) {
            out[i] = gradientD[i * d + i];
        }

        for (int b = 0; b < blockParameter.getNumBlocks(); ++b) {
            if (blockParameter.getBlockSizes()[b] != 2) {
                continue;
            }
            final int start = blockParameter.getBlockStarts()[b];
            out[upperBase + blockIndex] = gradientD[start * d + start + 1];
            out[lowerBase + blockIndex] = gradientD[(start + 1) * d + start];
            blockIndex++;
        }
    }

    private static void multiplyFlat(final double[] left,
                                     final double[] right,
                                     final double[] out,
                                     final int d) {
        for (int i = 0; i < d; ++i) {
            final int rowOffset = i * d;
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += left[rowOffset + k] * right[k * d + j];
                }
                out[rowOffset + j] = sum;
            }
        }
    }

    private static void multiplyTransposeLeftFlat(final double[] left,
                                                  final double[] right,
                                                  final double[] out,
                                                  final int d) {
        for (int i = 0; i < d; ++i) {
            final int rowOffset = i * d;
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += left[k * d + i] * right[k * d + j];
                }
                out[rowOffset + j] = sum;
            }
        }
    }

    private static void multiplyRightTransposeFlat(final double[] left,
                                                   final double[] right,
                                                   final double[] out,
                                                   final int d) {
        for (int i = 0; i < d; ++i) {
            final int rowOffset = i * d;
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                final int rightRowOffset = j * d;
                for (int k = 0; k < d; ++k) {
                    sum += left[rowOffset + k] * right[rightRowOffset + k];
                }
                out[rowOffset + j] = sum;
            }
        }
    }

}
