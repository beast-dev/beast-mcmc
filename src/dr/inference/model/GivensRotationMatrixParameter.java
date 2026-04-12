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
        final double[] work = new double[dim * dim];

        for (int k = 0; k < angleParameter.getDimension(); ++k) {
            final int i = pairI[k];
            final int j = pairJ[k];
            final double theta = angleParameter.getParameterValue(k);
            final double c = Math.cos(theta);
            final double s = Math.sin(theta);

            for (int row = 0; row < dim; ++row) {
                final int base = row * dim;
                for (int col = 0; col < dim; ++col) {
                    work[base + col] = matrixData[base + col];
                }
                work[base + i] = matrixData[base + i] * c + matrixData[base + j] * s;
                work[base + j] = -matrixData[base + i] * s + matrixData[base + j] * c;
            }

            System.arraycopy(work, 0, matrixData, 0, matrixData.length);
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
    public double[] pullBackGradient(final double[][] gradientWrtMatrix) {
        getParameterValue(0, 0);
        final int angleCount = angleParameter.getDimension();
        final double[] gradients = new double[angleCount];
        final double[][][] prefixes = new double[angleCount + 1][dim][dim];
        final double[][][] suffixes = new double[angleCount + 1][dim][dim];

        fillIdentity(prefixes[0], dim);
        for (int k = 0; k < angleCount; ++k) {
            final double[][] givens = buildGivens(k);
            multiply(prefixes[k], givens, prefixes[k + 1], dim);
        }

        fillIdentity(suffixes[angleCount], dim);
        for (int k = angleCount - 1; k >= 0; --k) {
            final double[][] givens = buildGivens(k);
            multiply(givens, suffixes[k + 1], suffixes[k], dim);
        }

        final double[][] tempLeft = new double[dim][dim];
        final double[][] temp = new double[dim][dim];
        for (int k = 0; k < angleCount; ++k) {
            final double[][] dG = buildGivensDerivative(k);
            transpose(prefixes[k], tempLeft, dim);
            multiply(tempLeft, gradientWrtMatrix, temp, dim);
            final double[][] suffixTranspose = new double[dim][dim];
            transpose(suffixes[k + 1], suffixTranspose, dim);
            multiply(temp, suffixTranspose, tempLeft, dim);

            double grad = 0.0;
            for (int i = 0; i < dim; ++i) {
                for (int j = 0; j < dim; ++j) {
                    grad += tempLeft[i][j] * dG[i][j];
                }
            }
            gradients[k] = grad;
        }

        return gradients;
    }

    public Parameter getAngleParameter() {
        return angleParameter;
    }

    @Override
    public boolean isConstrainedSymmetric() {
        return false;
    }

    private double[][] buildGivens(final int index) {
        final double[][] givens = new double[dim][dim];
        fillIdentity(givens, dim);
        final int i = pairI[index];
        final int j = pairJ[index];
        final double theta = angleParameter.getParameterValue(index);
        final double c = Math.cos(theta);
        final double s = Math.sin(theta);
        givens[i][i] = c;
        givens[i][j] = -s;
        givens[j][i] = s;
        givens[j][j] = c;
        return givens;
    }

    private double[][] buildGivensDerivative(final int index) {
        final double[][] deriv = new double[dim][dim];
        final int i = pairI[index];
        final int j = pairJ[index];
        final double theta = angleParameter.getParameterValue(index);
        final double c = Math.cos(theta);
        final double s = Math.sin(theta);
        deriv[i][i] = -s;
        deriv[i][j] = -c;
        deriv[j][i] = c;
        deriv[j][j] = -s;
        return deriv;
    }

    private static void fillIdentity(final double[] matrix, final int dimension) {
        for (int i = 0; i < dimension * dimension; ++i) {
            matrix[i] = 0.0;
        }
        for (int i = 0; i < dimension; ++i) {
            matrix[i * dimension + i] = 1.0;
        }
    }

    private static void fillIdentity(final double[][] matrix, final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                matrix[i][j] = i == j ? 1.0 : 0.0;
            }
        }
    }

    private static void transpose(final double[] in, final double[] out, final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                out[i * dimension + j] = in[j * dimension + i];
            }
        }
    }

    private static void transpose(final double[][] in, final double[][] out, final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                out[i][j] = in[j][i];
            }
        }
    }

    private static void multiply(final double[][] left,
                                 final double[][] right,
                                 final double[][] out,
                                 final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += left[i][k] * right[k][j];
                }
                out[i][j] = sum;
            }
        }
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
