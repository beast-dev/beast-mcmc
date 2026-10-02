package dr.inference.model;

import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

import java.util.Arrays;
import java.util.List;

/**
 * Invertible matrix parameterized as R = Q U.
 *
 * Q is supplied by an existing orthogonal matrix parametrization. U is unit
 * upper triangular, with free strict-upper entries.
 */
public final class QRUnitUpperMatrixParametrization extends InvertibleMatrixParametrization {

    public static final String NAME = "qrUnitUpperMatrixParametrization";

    private static final String ORTHOGONAL_MATRIX = "orthogonalMatrix";
    private static final String UNIT_UPPER = "unitUpper";

    private final MatrixParameter orthogonalMatrix;
    private final OrthogonalMatrixProvider orthogonalProvider;
    private final Parameter upperParameter;
    private final int[] upperRowOffsets;
    private final double[] qData;
    private final double[] qTransposeData;
    private final double[] unitUpperData;
    private final double[] unitUpperInverseData;
    private final double[] matrixData;
    private final double[] inverseData;
    private final double[] gradientQData;
    private final double[] gradientUpperData;

    public QRUnitUpperMatrixParametrization(final String name,
                                           final MatrixParameter orthogonalMatrix,
                                           final Parameter upperParameter) {
        super(name,
                requireSquareDimension(orthogonalMatrix),
                parameters(requireSquareDimension(orthogonalMatrix),
                        requireOrthogonalProvider(orthogonalMatrix),
                        upperParameter));
        this.orthogonalMatrix = orthogonalMatrix;
        this.orthogonalProvider = (OrthogonalMatrixProvider) orthogonalMatrix;
        this.upperParameter = upperParameter;
        this.upperRowOffsets = new int[dim];
        this.qData = new double[dim * dim];
        this.qTransposeData = new double[dim * dim];
        this.unitUpperData = new double[dim * dim];
        this.unitUpperInverseData = new double[dim * dim];
        this.matrixData = new double[dim * dim];
        this.inverseData = new double[dim * dim];
        this.gradientQData = new double[dim * dim];
        this.gradientUpperData = new double[getStrictUpperDimension(dim)];

        for (int row = 0; row < dim; row++) {
            upperRowOffsets[row] = StrictTriangularMatrixIndexing.upperRowOffset(dim, row);
        }
    }

    public MatrixParameter getOrthogonalMatrixParameter() {
        return orthogonalMatrix;
    }

    public OrthogonalMatrixProvider getOrthogonalProvider() {
        return orthogonalProvider;
    }

    public Parameter getOrthogonalParameter() {
        return orthogonalProvider.getOrthogonalParameter();
    }

    public Parameter getUpperParameter() {
        return upperParameter;
    }

    public static int getStrictUpperDimension(final int dimension) {
        return StrictTriangularMatrixIndexing.getStrictTriangularDimension(dimension);
    }

    public static int upperIndex(final int dimension, final int row, final int col) {
        return StrictTriangularMatrixIndexing.upperIndex(dimension, row, col);
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
        orthogonalProvider.fillOrthogonalMatrix(qData);
        orthogonalProvider.fillOrthogonalTranspose(qTransposeData);
        fillUnitUpperData();
        fillMatrixData();
        fillUnitUpperInverseData();
        fillInverseData();
    }

    @Override
    public double getLogAbsDeterminant() {
        return 0.0;
    }

    @Override
    public double getDeterminantSign() {
        return 1.0;
    }

    @Override
    public double getDeterminant() {
        return 1.0;
    }

    @Override
    public synchronized void fillInverse(final double[] inverseRowMajor) {
        if (inverseRowMajor == null || inverseRowMajor.length < dim * dim) {
            throw new IllegalArgumentException("inverseRowMajor must have length at least " + (dim * dim));
        }
        getParameterValue(0, 0);
        System.arraycopy(inverseData, 0, inverseRowMajor, 0, dim * dim);
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

        getParameterValue(0, 0);
        fillOrthogonalPullBackSource(gradientWrtMatrixRowMajor);
        fillUnitUpperPullBackSource(gradientWrtMatrixRowMajor);

        final int angleDimension = orthogonalProvider.getOrthogonalParameter().getDimension();
        orthogonalProvider.fillPullBackGradientFlat(gradientQData, dim, out, offset);
        System.arraycopy(gradientUpperData, 0, out, offset + angleDimension, gradientUpperData.length);
    }

    @Override
    public boolean isConstrainedSymmetric() {
        return false;
    }

    private void fillUnitUpperData() {
        Arrays.fill(unitUpperData, 0.0);
        for (int row = 0; row < dim; row++) {
            unitUpperData[row * dim + row] = 1.0;
            int upperIndex = upperRowOffsets[row];
            for (int col = row + 1; col < dim; col++) {
                unitUpperData[row * dim + col] = upperParameter.getParameterValue(upperIndex++);
            }
        }
    }

    private void fillMatrixData() {
        for (int row = 0; row < dim; row++) {
            final int rowOffset = row * dim;
            for (int col = 0; col < dim; col++) {
                double sum = qData[rowOffset + col];
                for (int k = 0; k < col; k++) {
                    sum += qData[rowOffset + k] * unitUpperData[k * dim + col];
                }
                matrixData[rowOffset + col] = sum;
            }
        }
    }

    private void fillUnitUpperInverseData() {
        Arrays.fill(unitUpperInverseData, 0.0);
        for (int row = dim - 1; row >= 0; row--) {
            unitUpperInverseData[row * dim + row] = 1.0;
            for (int col = row + 1; col < dim; col++) {
                double sum = 0.0;
                for (int k = row + 1; k <= col; k++) {
                    sum += unitUpperData[row * dim + k] * unitUpperInverseData[k * dim + col];
                }
                unitUpperInverseData[row * dim + col] = -sum;
            }
        }
    }

    private void fillInverseData() {
        for (int row = 0; row < dim; row++) {
            final int rowOffset = row * dim;
            for (int col = 0; col < dim; col++) {
                double sum = 0.0;
                for (int k = row; k < dim; k++) {
                    sum += unitUpperInverseData[rowOffset + k] * qTransposeData[k * dim + col];
                }
                inverseData[rowOffset + col] = sum;
            }
        }
    }

    private void fillOrthogonalPullBackSource(final double[] gradientR) {
        for (int row = 0; row < dim; row++) {
            final int rowOffset = row * dim;
            for (int col = 0; col < dim; col++) {
                double sum = 0.0;
                for (int k = col; k < dim; k++) {
                    sum += gradientR[rowOffset + k] * unitUpperData[col * dim + k];
                }
                gradientQData[rowOffset + col] = sum;
            }
        }
    }

    private void fillUnitUpperPullBackSource(final double[] gradientR) {
        int index = 0;
        for (int row = 0; row < dim - 1; row++) {
            for (int col = row + 1; col < dim; col++) {
                double sum = 0.0;
                for (int k = 0; k < dim; k++) {
                    sum += qData[k * dim + row] * gradientR[k * dim + col];
                }
                gradientUpperData[index++] = sum;
            }
        }
    }

    private static int requireSquareDimension(final MatrixParameter orthogonalMatrix) {
        if (orthogonalMatrix == null) {
            throw new IllegalArgumentException(ORTHOGONAL_MATRIX + " is required");
        }
        final int rows = orthogonalMatrix.getRowDimension();
        final int cols = orthogonalMatrix.getColumnDimension();
        if (rows != cols) {
            throw new IllegalArgumentException(ORTHOGONAL_MATRIX + " must be square");
        }
        return rows;
    }

    private static OrthogonalMatrixProvider requireOrthogonalProvider(final MatrixParameter orthogonalMatrix) {
        if (!(orthogonalMatrix instanceof OrthogonalMatrixProvider)) {
            throw new IllegalArgumentException(
                    ORTHOGONAL_MATRIX + " must implement " + OrthogonalMatrixProvider.class.getSimpleName());
        }
        return (OrthogonalMatrixProvider) orthogonalMatrix;
    }

    private static List<Parameter> parameters(final int dimension,
                                              final OrthogonalMatrixProvider orthogonalProvider,
                                              final Parameter upperParameter) {
        if (orthogonalProvider.getOrthogonalParameter() == null) {
            throw new IllegalArgumentException(ORTHOGONAL_MATRIX + " native parameter is required");
        }
        requireDimension(upperParameter, getStrictUpperDimension(dimension), UNIT_UPPER);
        return Arrays.asList(orthogonalProvider.getOrthogonalParameter(), upperParameter);
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

    public static final XMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public String getParserName() {
            return NAME;
        }

        @Override
        public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
            final String name = xo.hasId() ? xo.getId() : null;
            final MatrixParameter orthogonal = (MatrixParameter) xo.getElementFirstChild(ORTHOGONAL_MATRIX);
            final Parameter unitUpper = (Parameter) xo.getElementFirstChild(UNIT_UPPER);
            return new QRUnitUpperMatrixParametrization(name, orthogonal, unitUpper);
        }

        @Override
        public String getParserDescription() {
            return "An invertible matrix parameterized as Q U with unit upper triangular U.";
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return RULES;
        }

        @Override
        public Class getReturnType() {
            return QRUnitUpperMatrixParametrization.class;
        }
    };

    private static final XMLSyntaxRule[] RULES = new XMLSyntaxRule[]{
            new ElementRule(ORTHOGONAL_MATRIX, new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class)}),
            new ElementRule(UNIT_UPPER, new XMLSyntaxRule[]{new ElementRule(Parameter.class)})
    };
}
