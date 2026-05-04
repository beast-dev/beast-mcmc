package dr.inference.model;

import dr.xml.*;

import java.util.Collections;

/**
 * Orthogonal matrix parameterized by a sequence of Givens rotation angles.
 */
public final class GivensRotationMatrixParameter extends AbstractComputedCompoundMatrix
        implements MatrixParameterInterface, OrthogonalMatrixProvider {

    public static final String NAME = "givensRotationMatrixParameter";

    private static final String DIMENSION = "dimension";
    private static final String ANGLES = "angles";

    private final Parameter angleParameter;
    private final int[] pairI;
    private final int[] pairJ;
    private final double[] matrixData;
    private final double[] transposeData;
    private final double[] gradientMatrixScratch;
    private final double[] prefixes;
    private final double[] suffixes;
    private final double[] givensScratch;
    private final double[] givensDerivativeScratch;
    private final double[] tempLeftScratch;
    private final double[] tempScratch;
    private final double[] suffixTransposeScratch;

    public GivensRotationMatrixParameter(final String name,
                                         final int dimension,
                                         final Parameter angleParameter) {
        super(name, dimension, Collections.singletonList(angleParameter));
        final int expectedAngles = dimension * (dimension - 1) / 2;
        if (angleParameter.getDimension() != expectedAngles) {
            throw new IllegalArgumentException(
                    "Expected " + expectedAngles + " Givens angles for dimension " + dimension
                            + " but found " + angleParameter.getDimension());
        }
        this.angleParameter = angleParameter;
        this.pairI = new int[expectedAngles];
        this.pairJ = new int[expectedAngles];
        this.matrixData = new double[dimension * dimension];
        this.transposeData = new double[dimension * dimension];
        this.gradientMatrixScratch = new double[dimension * dimension];
        this.prefixes = new double[(expectedAngles + 1) * dimension * dimension];
        this.suffixes = new double[(expectedAngles + 1) * dimension * dimension];
        this.givensScratch = new double[dimension * dimension];
        this.givensDerivativeScratch = new double[dimension * dimension];
        this.tempLeftScratch = new double[dimension * dimension];
        this.tempScratch = new double[dimension * dimension];
        this.suffixTransposeScratch = new double[dimension * dimension];

        int index = 0;
        for (int i = 0; i < dimension - 1; ++i) {
            for (int j = i + 1; j < dimension; ++j) {
                pairI[index] = i;
                pairJ[index] = j;
                index++;
            }
        }
    }

    @Override
    protected double computeEntry(final int row, final int col) {
        return matrixData[row * dim + col];
    }

    @Override
    protected double getCachedValue(final int row, final int col) {
        return matrixData[row * dim + col];
    }

    @Override
    protected void updateCache() {
        fillIdentity(matrixData, dim);

        for (int k = 0; k < angleParameter.getDimension(); ++k) {
            final int i = pairI[k];
            final int j = pairJ[k];
            final double theta = angleParameter.getParameterValue(k);
            final double c = Math.cos(theta);
            final double s = Math.sin(theta);

            for (int row = 0; row < dim; ++row) {
                final int base = row * dim;
                final double oldI = matrixData[base + i];
                final double oldJ = matrixData[base + j];
                matrixData[base + i] = oldI * c + oldJ * s;
                matrixData[base + j] = -oldI * s + oldJ * c;
            }
        }

        transpose(matrixData, transposeData, dim);
    }

    @Override
    public Parameter getOrthogonalParameter() {
        return angleParameter;
    }

    @Override
    public void fillOrthogonalMatrix(final double[] outRowMajor) {
        getParameterValue(0, 0);
        System.arraycopy(matrixData, 0, outRowMajor, 0, matrixData.length);
    }

    @Override
    public void fillOrthogonalTranspose(final double[] outRowMajor) {
        getParameterValue(0, 0);
        System.arraycopy(transposeData, 0, outRowMajor, 0, transposeData.length);
    }

    @Override
    public void applyOrthogonal(final double[] in, final double[] out) {
        System.arraycopy(in, 0, out, 0, dim);
        for (int k = angleParameter.getDimension() - 1; k >= 0; --k) {
            final int i = pairI[k];
            final int j = pairJ[k];
            final double theta = angleParameter.getParameterValue(k);
            final double c = Math.cos(theta);
            final double s = Math.sin(theta);
            final double xi = out[i];
            final double xj = out[j];
            out[i] = c * xi - s * xj;
            out[j] = s * xi + c * xj;
        }
    }

    @Override
    public void applyOrthogonalTranspose(final double[] in, final double[] out) {
        System.arraycopy(in, 0, out, 0, dim);
        for (int k = 0; k < angleParameter.getDimension(); ++k) {
            final int i = pairI[k];
            final int j = pairJ[k];
            final double theta = angleParameter.getParameterValue(k);
            final double c = Math.cos(theta);
            final double s = Math.sin(theta);
            final double xi = out[i];
            final double xj = out[j];
            out[i] = c * xi + s * xj;
            out[j] = -s * xi + c * xj;
        }
    }

    @Override
    public synchronized double[] pullBackGradient(final double[][] gradientWrtMatrix) {
        final int angleCount = angleParameter.getDimension();
        final double[] gradients = new double[angleCount];
        for (int i = 0; i < dim; ++i) {
            if (gradientWrtMatrix[i] == null || gradientWrtMatrix[i].length != dim) {
                throw new IllegalArgumentException("gradientWrtMatrix must be square with dimension " + dim);
            }
            System.arraycopy(gradientWrtMatrix[i], 0, gradientMatrixScratch, i * dim, dim);
        }
        fillPullBackGradientFlat(gradientMatrixScratch, dim, gradients);
        return gradients;
    }

    @Override
    public synchronized void fillPullBackGradientFlat(final double[] gradientWrtMatrixRowMajor,
                                         final int dimension,
                                         final double[] out,
                                         final int offset) {
        if (dimension != dim || gradientWrtMatrixRowMajor == null
                || gradientWrtMatrixRowMajor.length != dim * dim) {
            throw new IllegalArgumentException(
                    "gradientWrtMatrixRowMajor must have length " + (dim * dim));
        }
        final int angleCount = angleParameter.getDimension();
        if (out == null || out.length < offset + angleCount) {
            throw new IllegalArgumentException("out must have room for " + angleCount + " entries");
        }
        if (!Boolean.getBoolean("beast.experimental.projectedFlatGivensPullback")) {
            fillPullBackGradientFlatExact(gradientWrtMatrixRowMajor, out, offset);
        } else {
            fillPullBackGradientFlatProjected(gradientWrtMatrixRowMajor, out, offset);
        }
    }

    private void fillPullBackGradientFlatProjected(final double[] gradientWrtMatrixRowMajor,
                                                   final double[] out,
                                                   final int offset) {
        getParameterValue(0, 0);
        final int angleCount = angleParameter.getDimension();
        final int matrixLength = dim * dim;

        fillIdentity(prefixes, 0, dim);
        for (int k = 0; k < angleCount; ++k) {
            fillGivens(k, givensScratch);
            multiplyFlat(prefixes, k * matrixLength, givensScratch, 0,
                    prefixes, (k + 1) * matrixLength, dim);
        }

        fillIdentity(suffixes, angleCount * matrixLength, dim);
        for (int k = angleCount - 1; k >= 0; --k) {
            fillGivens(k, givensScratch);
            multiplyFlat(givensScratch, 0, suffixes, (k + 1) * matrixLength,
                    suffixes, k * matrixLength, dim);
        }

        for (int k = 0; k < angleCount; ++k) {
            final int i = pairI[k];
            final int j = pairJ[k];
            final double theta = angleParameter.getParameterValue(k);
            final double c = Math.cos(theta);
            final double s = Math.sin(theta);
            final int prefixOffset = k * matrixLength;
            final int suffixOffset = (k + 1) * matrixLength;

            out[offset + k] =
                    -s * projectedGradientEntry(prefixes, prefixOffset, gradientWrtMatrixRowMajor,
                            suffixes, suffixOffset, i, i, dim)
                            - c * projectedGradientEntry(prefixes, prefixOffset, gradientWrtMatrixRowMajor,
                            suffixes, suffixOffset, i, j, dim)
                            + c * projectedGradientEntry(prefixes, prefixOffset, gradientWrtMatrixRowMajor,
                            suffixes, suffixOffset, j, i, dim)
                            - s * projectedGradientEntry(prefixes, prefixOffset, gradientWrtMatrixRowMajor,
                            suffixes, suffixOffset, j, j, dim);
        }
    }

    private void fillPullBackGradientFlatExact(final double[] gradientWrtMatrixRowMajor,
                                               final double[] out,
                                               final int offset) {
        getParameterValue(0, 0);
        final int angleCount = angleParameter.getDimension();
        final int matrixLength = dim * dim;

        fillIdentity(prefixes, 0, dim);
        for (int k = 0; k < angleCount; ++k) {
            fillGivens(k, givensScratch);
            multiplyFlat(prefixes, k * matrixLength, givensScratch, 0,
                    prefixes, (k + 1) * matrixLength, dim);
        }

        fillIdentity(suffixes, angleCount * matrixLength, dim);
        for (int k = angleCount - 1; k >= 0; --k) {
            fillGivens(k, givensScratch);
            multiplyFlat(givensScratch, 0, suffixes, (k + 1) * matrixLength,
                    suffixes, k * matrixLength, dim);
        }

        for (int k = 0; k < angleCount; ++k) {
            fillGivensDerivative(k, givensDerivativeScratch);
            transposeFlat(prefixes, k * matrixLength, tempLeftScratch, 0, dim);
            multiplyFlat(tempLeftScratch, 0, gradientWrtMatrixRowMajor, 0, tempScratch, 0, dim);
            transposeFlat(suffixes, (k + 1) * matrixLength, suffixTransposeScratch, 0, dim);
            multiplyFlat(tempScratch, 0, suffixTransposeScratch, 0, tempLeftScratch, 0, dim);

            double grad = 0.0;
            for (int i = 0; i < matrixLength; ++i) {
                grad += tempLeftScratch[i] * givensDerivativeScratch[i];
            }
            out[offset + k] = grad;
        }
    }

    public Parameter getAngleParameter() {
        return angleParameter;
    }

    @Override
    public boolean isConstrainedSymmetric() {
        return false;
    }

    private static void fillIdentity(final double[] matrix, final int dimension) {
        for (int i = 0; i < dimension * dimension; ++i) {
            matrix[i] = 0.0;
        }
        for (int i = 0; i < dimension; ++i) {
            matrix[i * dimension + i] = 1.0;
        }
    }

    private static void fillIdentity(final double[] matrix,
                                     final int offset,
                                     final int dimension) {
        for (int i = 0; i < dimension * dimension; ++i) {
            matrix[offset + i] = 0.0;
        }
        for (int i = 0; i < dimension; ++i) {
            matrix[offset + i * dimension + i] = 1.0;
        }
    }

    private void fillGivens(final int index, final double[] out) {
        fillIdentity(out, dim);
        final int i = pairI[index];
        final int j = pairJ[index];
        final double theta = angleParameter.getParameterValue(index);
        final double c = Math.cos(theta);
        final double s = Math.sin(theta);
        out[i * dim + i] = c;
        out[i * dim + j] = -s;
        out[j * dim + i] = s;
        out[j * dim + j] = c;
    }

    private void fillGivensDerivative(final int index, final double[] out) {
        for (int i = 0; i < out.length; ++i) {
            out[i] = 0.0;
        }
        final int i = pairI[index];
        final int j = pairJ[index];
        final double theta = angleParameter.getParameterValue(index);
        final double c = Math.cos(theta);
        final double s = Math.sin(theta);
        out[i * dim + i] = -s;
        out[i * dim + j] = -c;
        out[j * dim + i] = c;
        out[j * dim + j] = -s;
    }

    private static void transpose(final double[] in, final double[] out, final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                out[i * dimension + j] = in[j * dimension + i];
            }
        }
    }

    private static void transposeFlat(final double[] in,
                                      final int inOffset,
                                      final double[] out,
                                      final int outOffset,
                                      final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                out[outOffset + i * dimension + j] = in[inOffset + j * dimension + i];
            }
        }
    }

    private static void multiplyFlat(final double[] left,
                                     final int leftOffset,
                                     final double[] right,
                                     final int rightOffset,
                                     final double[] out,
                                     final int outOffset,
                                     final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += left[leftOffset + rowOffset + k]
                            * right[rightOffset + k * dimension + j];
                }
                out[outOffset + rowOffset + j] = sum;
            }
        }
    }

    private static double projectedGradientEntry(final double[] prefix,
                                                 final int prefixOffset,
                                                 final double[] gradient,
                                                 final double[] suffix,
                                                 final int suffixOffset,
                                                 final int row,
                                                 final int col,
                                                 final int dimension) {
        double sum = 0.0;
        for (int p = 0; p < dimension; ++p) {
            final double prefixValue = prefix[prefixOffset + p * dimension + row];
            final int gradientRowOffset = p * dimension;
            for (int q = 0; q < dimension; ++q) {
                sum += prefixValue * gradient[gradientRowOffset + q]
                        * suffix[suffixOffset + col * dimension + q];
            }
        }
        return sum;
    }

    public static final XMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public String getParserName() {
            return NAME;
        }

        @Override
        public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
            final String name = xo.hasId() ? xo.getId() : null;
            final int dimension = xo.getIntegerAttribute(DIMENSION);
            final Parameter angles = (Parameter) xo.getElementFirstChild(ANGLES);
            return new GivensRotationMatrixParameter(name, dimension, angles);
        }

        @Override
        public String getParserDescription() {
            return "An orthogonal matrix parameterized by Givens angles.";
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return RULES;
        }

        @Override
        public Class getReturnType() {
            return GivensRotationMatrixParameter.class;
        }
    };

    private static final XMLSyntaxRule[] RULES = new XMLSyntaxRule[]{
            AttributeRule.newIntegerRule(DIMENSION),
            new ElementRule(ANGLES, new XMLSyntaxRule[]{new ElementRule(Parameter.class)})
    };
}
