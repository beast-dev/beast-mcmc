package dr.inference.model;

import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

import java.util.Arrays;
import java.util.List;

/**
 * Invertible matrix parameterized as R = L D U.
 *
 * L and U are unit lower/upper triangular matrices with unconstrained strict
 * triangular entries. D has fixed signs and exp(logDiagonal) magnitudes, so no
 * proposed state can make the matrix singular.
 */
public final class LDUInvertibleMatrixParametrization extends InvertibleMatrixParametrization {

    public static final String NAME = "lduInvertibleMatrixParametrization";

    private static final String DIMENSION = "dimension";
    private static final String LOWER = "lower";
    private static final String LOG_DIAGONAL = "logDiagonal";
    private static final String UPPER = "upper";
    private static final String DIAGONAL_SIGNS = "diagonalSigns";

    private final Parameter lowerParameter;
    private final Parameter logDiagonalParameter;
    private final Parameter upperParameter;
    private final double[] diagonalSigns;
    private final int[] lowerRowOffsets;
    private final int[] upperRowOffsets;
    private final double[] matrixData;
    private final double[] diagonalData;
    private final double[] gradientLowerData;
    private final double[] gradientLogDiagonalData;
    private final double[] gradientUpperData;
    private final double[] inverseYScratch;
    private final double[] inverseZScratch;
    private final double[] inverseXScratch;

    public LDUInvertibleMatrixParametrization(final String name,
                                             final int dimension,
                                             final Parameter lowerParameter,
                                             final Parameter logDiagonalParameter,
                                             final Parameter upperParameter) {
        this(name, dimension, lowerParameter, logDiagonalParameter, upperParameter, null);
    }

    public LDUInvertibleMatrixParametrization(final String name,
                                             final int dimension,
                                             final Parameter lowerParameter,
                                             final Parameter logDiagonalParameter,
                                             final Parameter upperParameter,
                                             final double[] diagonalSigns) {
        super(name, dimension, parameters(dimension, lowerParameter, logDiagonalParameter, upperParameter));
        this.lowerParameter = lowerParameter;
        this.logDiagonalParameter = logDiagonalParameter;
        this.upperParameter = upperParameter;
        this.diagonalSigns = normalizeSigns(dimension, diagonalSigns);
        this.lowerRowOffsets = new int[dimension];
        this.upperRowOffsets = new int[dimension];
        this.matrixData = new double[dimension * dimension];
        this.diagonalData = new double[dimension];
        this.gradientLowerData = new double[getStrictTriangularDimension(dimension)];
        this.gradientLogDiagonalData = new double[dimension];
        this.gradientUpperData = new double[getStrictTriangularDimension(dimension)];
        this.inverseYScratch = new double[dimension];
        this.inverseZScratch = new double[dimension];
        this.inverseXScratch = new double[dimension];

        for (int row = 0; row < dimension; row++) {
            lowerRowOffsets[row] = row * (row - 1) / 2;
            upperRowOffsets[row] = row * (2 * dimension - row - 1) / 2;
        }
    }

    public Parameter getLowerParameter() {
        return lowerParameter;
    }

    public Parameter getLogDiagonalParameter() {
        return logDiagonalParameter;
    }

    public Parameter getUpperParameter() {
        return upperParameter;
    }

    public double getDiagonalSign(final int index) {
        return diagonalSigns[index];
    }

    public static int getStrictTriangularDimension(final int dimension) {
        return dimension * (dimension - 1) / 2;
    }

    public static int lowerIndex(final int row, final int col) {
        if (row <= col) {
            throw new IllegalArgumentException("lowerIndex requires row > col");
        }
        return row * (row - 1) / 2 + col;
    }

    public static int upperIndex(final int dimension, final int row, final int col) {
        if (row >= col) {
            throw new IllegalArgumentException("upperIndex requires row < col");
        }
        return row * (2 * dimension - row - 1) / 2 + (col - row - 1);
    }

    @Override
    protected double computeEntry(final int row, final int col) {
        fillDiagonalData();
        return computeEntryWithCachedDiagonal(row, col);
    }

    private double computeEntryWithCachedDiagonal(final int row, final int col) {
        double sum = 0.0;
        final int maxSharedIndex = Math.min(row, col);
        final int lowerRowOffset = lowerRowOffsets[row];
        for (int k = 0; k <= maxSharedIndex; k++) {
            final double lower = (row == k)
                    ? 1.0
                    : lowerParameter.getParameterValue(lowerRowOffset + k);
            final double upper = (k == col)
                    ? 1.0
                    : upperParameter.getParameterValue(upperRowOffsets[k] + col - k - 1);
            sum += lower * diagonalData[k] * upper;
        }
        return sum;
    }

    @Override
    protected double getCachedValue(final int row, final int col) {
        return matrixData[row * dim + col];
    }

    @Override
    protected void updateCache() {
        fillDiagonalData();
        for (int row = 0; row < dim; row++) {
            for (int col = 0; col < dim; col++) {
                matrixData[row * dim + col] = computeEntryWithCachedDiagonal(row, col);
            }
        }
    }

    @Override
    public double getLogAbsDeterminant() {
        double logAbsDeterminant = 0.0;
        for (int i = 0; i < dim; i++) {
            logAbsDeterminant += logDiagonalParameter.getParameterValue(i);
        }
        return logAbsDeterminant;
    }

    @Override
    public double getDeterminantSign() {
        double sign = 1.0;
        for (double diagonalSign : diagonalSigns) {
            sign *= diagonalSign;
        }
        return sign;
    }

    @Override
    public synchronized void fillInverse(final double[] inverseRowMajor) {
        if (inverseRowMajor == null || inverseRowMajor.length < dim * dim) {
            throw new IllegalArgumentException("inverseRowMajor must have length at least " + (dim * dim));
        }

        fillDiagonalData();

        for (int col = 0; col < dim; col++) {
            for (int i = 0; i < dim; i++) {
                double value = (i == col) ? 1.0 : 0.0;
                final int lowerRowOffset = lowerRowOffsets[i];
                for (int k = 0; k < i; k++) {
                    value -= lowerParameter.getParameterValue(lowerRowOffset + k) * inverseYScratch[k];
                }
                inverseYScratch[i] = value;
            }

            for (int i = 0; i < dim; i++) {
                inverseZScratch[i] = inverseYScratch[i] / diagonalData[i];
            }

            for (int i = dim - 1; i >= 0; i--) {
                double value = inverseZScratch[i];
                int upperIndex = upperRowOffsets[i];
                for (int k = i + 1; k < dim; k++) {
                    value -= upperParameter.getParameterValue(upperIndex++) * inverseXScratch[k];
                }
                inverseXScratch[i] = value;
            }

            for (int row = 0; row < dim; row++) {
                inverseRowMajor[row * dim + col] = inverseXScratch[row];
            }
        }
    }

    @Override
    public synchronized void fillPullBackGradientFlat(final double[] gradientWrtMatrixRowMajor,
                                                      final int dimension,
                                                      final double[] out,
                                                      final int offset) {
        if (dimension != dim
                || gradientWrtMatrixRowMajor == null
                || gradientWrtMatrixRowMajor.length != dim * dim) {
            throw new IllegalArgumentException(
                    "gradientWrtMatrixRowMajor must have length " + (dim * dim));
        }
        if (out == null || out.length < offset + getNativeDimension()) {
            throw new IllegalArgumentException(
                    "out must have room for " + getNativeDimension() + " entries at offset " + offset);
        }

        fillDiagonalData();
        fillLowerPullBackGradient(gradientWrtMatrixRowMajor);
        fillLogDiagonalPullBackGradient(gradientWrtMatrixRowMajor);
        fillUpperPullBackGradient(gradientWrtMatrixRowMajor);

        int outOffset = offset;
        System.arraycopy(gradientLowerData, 0, out, outOffset, gradientLowerData.length);
        outOffset += gradientLowerData.length;
        System.arraycopy(gradientLogDiagonalData, 0, out, outOffset, gradientLogDiagonalData.length);
        outOffset += gradientLogDiagonalData.length;
        System.arraycopy(gradientUpperData, 0, out, outOffset, gradientUpperData.length);
    }

    private double diagonalEntry(final int index) {
        return diagonalSigns[index] * Math.exp(logDiagonalParameter.getParameterValue(index));
    }

    private void fillDiagonalData() {
        for (int i = 0; i < dim; i++) {
            diagonalData[i] = diagonalEntry(i);
        }
    }

    private void fillLowerPullBackGradient(final double[] gradientWrtMatrixRowMajor) {
        int index = 0;
        for (int row = 1; row < dim; row++) {
            final int gradientRowOffset = row * dim;
            for (int col = 0; col < row; col++) {
                double rowSum = gradientWrtMatrixRowMajor[gradientRowOffset + col];
                int upperIndex = upperRowOffsets[col];
                for (int k = col + 1; k < dim; k++) {
                    rowSum += gradientWrtMatrixRowMajor[gradientRowOffset + k]
                            * upperParameter.getParameterValue(upperIndex++);
                }
                gradientLowerData[index++] = diagonalData[col] * rowSum;
            }
        }
    }

    private void fillLogDiagonalPullBackGradient(final double[] gradientWrtMatrixRowMajor) {
        for (int diag = 0; diag < dim; diag++) {
            final int upperRowOffset = upperRowOffsets[diag];
            int gradientRowOffset = diag * dim;
            double gradientD = gradientWrtMatrixRowMajor[gradientRowOffset + diag];
            int upperIndex = upperRowOffset;
            for (int col = diag + 1; col < dim; col++) {
                gradientD += gradientWrtMatrixRowMajor[gradientRowOffset + col]
                        * upperParameter.getParameterValue(upperIndex++);
            }
            for (int row = diag + 1; row < dim; row++) {
                final double lower = lowerParameter.getParameterValue(lowerRowOffsets[row] + diag);
                gradientRowOffset = row * dim;
                double rowSum = gradientWrtMatrixRowMajor[gradientRowOffset + diag];
                upperIndex = upperRowOffset;
                for (int col = diag + 1; col < dim; col++) {
                    rowSum += gradientWrtMatrixRowMajor[gradientRowOffset + col]
                            * upperParameter.getParameterValue(upperIndex++);
                }
                gradientD += lower * rowSum;
            }
            gradientLogDiagonalData[diag] = gradientD * diagonalData[diag];
        }
    }

    private void fillUpperPullBackGradient(final double[] gradientWrtMatrixRowMajor) {
        int index = 0;
        for (int row = 0; row < dim - 1; row++) {
            final int gradientRowOffset = row * dim;
            for (int col = row + 1; col < dim; col++) {
                double columnSum = gradientWrtMatrixRowMajor[gradientRowOffset + col];
                for (int k = row + 1; k < dim; k++) {
                    columnSum += lowerParameter.getParameterValue(lowerRowOffsets[k] + row)
                            * gradientWrtMatrixRowMajor[k * dim + col];
                }
                gradientUpperData[index++] = diagonalData[row] * columnSum;
            }
        }
    }

    private static List<Parameter> parameters(final int dimension,
                                              final Parameter lowerParameter,
                                              final Parameter logDiagonalParameter,
                                              final Parameter upperParameter) {
        final int triangularDimension = getStrictTriangularDimension(dimension);
        requireDimension(lowerParameter, triangularDimension, LOWER);
        requireDimension(logDiagonalParameter, dimension, LOG_DIAGONAL);
        requireDimension(upperParameter, triangularDimension, UPPER);
        return Arrays.asList(lowerParameter, logDiagonalParameter, upperParameter);
    }

    private static void requireDimension(final Parameter parameter,
                                         final int expectedDimension,
                                         final String name) {
        if (parameter == null) {
            throw new IllegalArgumentException(name + " parameter is required");
        }
        if (parameter.getDimension() != expectedDimension) {
            throw new IllegalArgumentException(name + " parameter must have dimension "
                    + expectedDimension + " but found " + parameter.getDimension());
        }
    }

    private static double[] normalizeSigns(final int dimension, final double[] diagonalSigns) {
        final double[] signs = new double[dimension];
        if (diagonalSigns == null) {
            Arrays.fill(signs, 1.0);
            return signs;
        }
        if (diagonalSigns.length != dimension) {
            throw new IllegalArgumentException(DIAGONAL_SIGNS + " must have length " + dimension);
        }
        for (int i = 0; i < dimension; i++) {
            if (diagonalSigns[i] == 0.0) {
                throw new IllegalArgumentException(DIAGONAL_SIGNS + " entries must be non-zero");
            }
            signs[i] = diagonalSigns[i] < 0.0 ? -1.0 : 1.0;
        }
        return signs;
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
            final Parameter lower = (Parameter) xo.getElementFirstChild(LOWER);
            final Parameter logDiagonal = (Parameter) xo.getElementFirstChild(LOG_DIAGONAL);
            final Parameter upper = (Parameter) xo.getElementFirstChild(UPPER);
            final double[] diagonalSigns = xo.hasAttribute(DIAGONAL_SIGNS)
                    ? xo.getDoubleArrayAttribute(DIAGONAL_SIGNS)
                    : null;

            return new LDUInvertibleMatrixParametrization(
                    name, dimension, lower, logDiagonal, upper, diagonalSigns);
        }

        @Override
        public String getParserDescription() {
            return "An invertible matrix parameterized as L D U with exp(logDiagonal) pivots.";
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return RULES;
        }

        @Override
        public Class getReturnType() {
            return LDUInvertibleMatrixParametrization.class;
        }
    };

    private static final XMLSyntaxRule[] RULES = new XMLSyntaxRule[]{
            AttributeRule.newIntegerRule(DIMENSION),
            AttributeRule.newDoubleArrayRule(DIAGONAL_SIGNS, true),
            new ElementRule(LOWER, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(LOG_DIAGONAL, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(UPPER, new XMLSyntaxRule[]{new ElementRule(Parameter.class)})
    };
}
