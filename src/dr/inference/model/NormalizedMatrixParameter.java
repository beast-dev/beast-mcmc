package dr.inference.model;

import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

import java.util.Collections;

/**
 * Computed dense matrix parameter that normalizes a raw square matrix either by
 * one global Frobenius norm or independently by columns.
 */
public class NormalizedMatrixParameter extends AbstractComputedCompoundMatrix
        implements NativeMatrixGradientProvider {

    public static final String NAME = "normalizedMatrixParameter";

    public enum Normalization {
        GLOBAL("global"),
        COLUMNS("columns");

        private final String xmlName;

        Normalization(final String xmlName) {
            this.xmlName = xmlName;
        }

        public String getXmlName() {
            return xmlName;
        }

        public static Normalization parse(final String value) {
            for (final Normalization normalization : values()) {
                if (normalization.xmlName.equalsIgnoreCase(value)) {
                    return normalization;
                }
            }
            throw new IllegalArgumentException(
                    "Unsupported normalization mode: " + value
                            + " (expected global or columns)");
        }
    }

    private static final String NORMALIZATION = "normalization";
    private static final String MIN_NORM = "minNorm";
    public static final double DEFAULT_MIN_NORM = 1.0e-12;

    private final MatrixParameter rawMatrix;
    private final Normalization normalization;
    private final double minNorm;
    private final double[] normalizedData;
    private final double[] columnNorms;
    private double globalNorm;

    public NormalizedMatrixParameter(final String name,
                                     final MatrixParameter rawMatrix,
                                     final Normalization normalization) {
        this(name, rawMatrix, normalization, DEFAULT_MIN_NORM);
    }

    public NormalizedMatrixParameter(final String name,
                                     final MatrixParameter rawMatrix,
                                     final Normalization normalization,
                                     final double minNorm) {
        super(name, requireSquareDimension(rawMatrix), Collections.<Parameter>singletonList(rawMatrix));
        if (normalization == null) {
            throw new IllegalArgumentException(NORMALIZATION + " is required");
        }
        if (minNorm <= 0.0 || Double.isNaN(minNorm)) {
            throw new IllegalArgumentException(MIN_NORM + " must be positive");
        }
        this.rawMatrix = rawMatrix;
        this.normalization = normalization;
        this.minNorm = minNorm;
        this.normalizedData = new double[dim * dim];
        this.columnNorms = new double[dim];
        this.globalNorm = Double.NaN;
    }

    public MatrixParameter getRawMatrixParameter() {
        return rawMatrix;
    }

    public Normalization getNormalization() {
        return normalization;
    }

    public double getMinNorm() {
        return minNorm;
    }

    public double getGlobalNorm() {
        getParameterValue(0, 0);
        return globalNorm;
    }

    public double getColumnNorm(final int column) {
        getParameterValue(0, 0);
        return columnNorms[column];
    }

    @Override
    protected double computeEntry(final int row, final int col) {
        final double norm = normalization == Normalization.GLOBAL
                ? computeGlobalNorm()
                : computeColumnNorm(col);
        return rawMatrix.getParameterValue(row, col) / norm;
    }

    @Override
    protected double getCachedValue(final int row, final int col) {
        return normalizedData[row * dim + col];
    }

    @Override
    protected void updateCache() {
        if (normalization == Normalization.GLOBAL) {
            globalNorm = computeGlobalNorm();
            for (int col = 0; col < dim; col++) {
                columnNorms[col] = computeColumnNormUnchecked(col);
                for (int row = 0; row < dim; row++) {
                    normalizedData[row * dim + col] = rawMatrix.getParameterValue(row, col) / globalNorm;
                }
            }
            return;
        }

        globalNorm = computeGlobalNormUnchecked();
        for (int col = 0; col < dim; col++) {
            final double norm = computeColumnNorm(col);
            columnNorms[col] = norm;
            for (int row = 0; row < dim; row++) {
                normalizedData[row * dim + col] = rawMatrix.getParameterValue(row, col) / norm;
            }
        }
    }

    @Override
    public Parameter getNativeParameter() {
        return rawMatrix;
    }

    @Override
    public boolean supportsNativeParameter(final Parameter parameter) {
        return parameter == rawMatrix;
    }

    @Override
    public int getPullBackGradientDimension(final Parameter parameter) {
        if (supportsNativeParameter(parameter)) {
            return rawMatrix.getDimension();
        }
        throw unsupportedNativeParameter(parameter);
    }

    @Override
    public double[] pullBackGradientForParameter(final Parameter parameter,
                                                 final double[] gradientWrtMatrixRowMajor,
                                                 final int dimension) {
        final double[] gradient = new double[getPullBackGradientDimension(parameter)];
        fillPullBackGradientForParameter(parameter, gradientWrtMatrixRowMajor, dimension, gradient, 0);
        return gradient;
    }

    @Override
    public synchronized void fillPullBackGradientForParameter(final Parameter parameter,
                                                              final double[] gradientWrtMatrixRowMajor,
                                                              final int dimension,
                                                              final double[] out,
                                                              final int offset) {
        final int outDimension = getPullBackGradientDimension(parameter);
        if (dimension != dim
                || gradientWrtMatrixRowMajor == null
                || gradientWrtMatrixRowMajor.length != dim * dim) {
            throw new IllegalArgumentException(
                    "gradientWrtMatrixRowMajor must have length " + (dim * dim));
        }
        if (out == null || out.length < offset + outDimension) {
            throw new IllegalArgumentException(
                    "out must have room for " + outDimension + " entries at offset " + offset);
        }

        getParameterValue(0, 0);
        if (normalization == Normalization.GLOBAL) {
            fillGlobalPullBack(gradientWrtMatrixRowMajor, out, offset);
        } else {
            fillColumnPullBack(gradientWrtMatrixRowMajor, out, offset);
        }
    }

    @Override
    public boolean isConstrainedSymmetric() {
        return false;
    }

    private void fillGlobalPullBack(final double[] gradientWrtMatrixRowMajor,
                                    final double[] out,
                                    final int offset) {
        double radialComponent = 0.0;
        for (int col = 0; col < dim; col++) {
            for (int row = 0; row < dim; row++) {
                final int matrixIndex = row * dim + col;
                radialComponent += normalizedData[matrixIndex] * gradientWrtMatrixRowMajor[matrixIndex];
            }
        }

        final double inverseNorm = 1.0 / globalNorm;
        for (int col = 0; col < dim; col++) {
            for (int row = 0; row < dim; row++) {
                final int matrixIndex = row * dim + col;
                final int nativeIndex = offset + col * dim + row;
                out[nativeIndex] = (gradientWrtMatrixRowMajor[matrixIndex]
                        - normalizedData[matrixIndex] * radialComponent) * inverseNorm;
            }
        }
    }

    private void fillColumnPullBack(final double[] gradientWrtMatrixRowMajor,
                                    final double[] out,
                                    final int offset) {
        for (int col = 0; col < dim; col++) {
            double radialComponent = 0.0;
            for (int row = 0; row < dim; row++) {
                final int matrixIndex = row * dim + col;
                radialComponent += normalizedData[matrixIndex] * gradientWrtMatrixRowMajor[matrixIndex];
            }

            final double inverseNorm = 1.0 / columnNorms[col];
            for (int row = 0; row < dim; row++) {
                final int matrixIndex = row * dim + col;
                final int nativeIndex = offset + col * dim + row;
                out[nativeIndex] = (gradientWrtMatrixRowMajor[matrixIndex]
                        - normalizedData[matrixIndex] * radialComponent) * inverseNorm;
            }
        }
    }

    private double computeGlobalNorm() {
        final double norm = computeGlobalNormUnchecked();
        if (norm <= minNorm || Double.isNaN(norm)) {
            throw new IllegalArgumentException(
                    rawMatrix.getId() + " has Frobenius norm " + norm
                            + ", which is not above " + minNorm);
        }
        return norm;
    }

    private double computeGlobalNormUnchecked() {
        double sumSquares = 0.0;
        for (int col = 0; col < dim; col++) {
            for (int row = 0; row < dim; row++) {
                final double value = rawMatrix.getParameterValue(row, col);
                sumSquares += value * value;
            }
        }
        return Math.sqrt(sumSquares);
    }

    private double computeColumnNorm(final int col) {
        final double norm = computeColumnNormUnchecked(col);
        if (norm <= minNorm || Double.isNaN(norm)) {
            throw new IllegalArgumentException(
                    "Column " + col + " of " + rawMatrix.getId()
                            + " has norm " + norm
                            + ", which is not above " + minNorm);
        }
        return norm;
    }

    private double computeColumnNormUnchecked(final int col) {
        double sumSquares = 0.0;
        for (int row = 0; row < dim; row++) {
            final double value = rawMatrix.getParameterValue(row, col);
            sumSquares += value * value;
        }
        return Math.sqrt(sumSquares);
    }

    private IllegalArgumentException unsupportedNativeParameter(final Parameter parameter) {
        return new IllegalArgumentException("Unsupported native parameter: "
                + (parameter == null ? "null" : parameter.getId()));
    }

    private static int requireSquareDimension(final MatrixParameter rawMatrix) {
        if (rawMatrix == null) {
            throw new IllegalArgumentException("raw matrix is required");
        }
        final int rows = rawMatrix.getRowDimension();
        final int cols = rawMatrix.getColumnDimension();
        if (rows != cols) {
            throw new IllegalArgumentException("raw matrix must be square");
        }
        return rows;
    }

    public static final XMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public String getParserName() {
            return NAME;
        }

        @Override
        public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
            final String name = xo.hasId() ? xo.getId() : null;
            final MatrixParameter raw = (MatrixParameter) xo.getChild(MatrixParameter.class);
            final String mode = xo.getStringAttribute(NORMALIZATION);
            final double minNorm = xo.getAttribute(MIN_NORM, DEFAULT_MIN_NORM);
            return new NormalizedMatrixParameter(
                    name,
                    raw,
                    Normalization.parse(mode),
                    minNorm);
        }

        @Override
        public String getParserDescription() {
            return "A dense matrix parameter normalized globally or by columns.";
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return RULES;
        }

        @Override
        public Class getReturnType() {
            return NormalizedMatrixParameter.class;
        }
    };

    private static final XMLSyntaxRule[] RULES = new XMLSyntaxRule[]{
            AttributeRule.newStringRule(NORMALIZATION, false),
            AttributeRule.newDoubleRule(MIN_NORM, true),
            new ElementRule(MatrixParameter.class)
    };
}
