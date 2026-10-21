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
 * Computed square matrix equal to one scalar parameter times the identity.
 */
public final class ScaledIdentityMatrixParameter extends AbstractComputedCompoundMatrix {

    public static final String NAME = "scaledIdentityMatrixParameter";
    private static final String DIMENSION = "dimension";

    private final Parameter scaleParameter;
    private final double[] values;

    public ScaledIdentityMatrixParameter(final String name,
                                         final int dimension,
                                         final Parameter scaleParameter) {
        super(name, requireDimension(dimension), Collections.<Parameter>singletonList(requireScalar(scaleParameter)));
        this.scaleParameter = scaleParameter;
        this.values = new double[dimension * dimension];
    }

    public Parameter getScaleParameter() {
        return scaleParameter;
    }

    @Override
    protected double computeEntry(final int row, final int col) {
        return row == col ? scaleParameter.getParameterValue(0) : 0.0;
    }

    @Override
    protected double getCachedValue(final int row, final int col) {
        return values[row * dim + col];
    }

    @Override
    protected void updateCache() {
        final double scale = scaleParameter.getParameterValue(0);
        for (int index = 0; index < values.length; ++index) {
            values[index] = 0.0;
        }
        for (int index = 0; index < dim; ++index) {
            values[index * dim + index] = scale;
        }
    }

    @Override
    public double[] getColumnValues(final int col) {
        getParameterValue(0, 0);
        final double[] column = new double[dim];
        for (int row = 0; row < dim; ++row) {
            column[row] = getCachedValue(row, col);
        }
        return column;
    }

    @Override
    public boolean isConstrainedSymmetric() {
        return true;
    }

    private static Parameter requireScalar(final Parameter parameter) {
        if (parameter == null) {
            throw new IllegalArgumentException("scale parameter is required");
        }
        if (parameter.getDimension() != 1) {
            throw new IllegalArgumentException(
                    "scaled identity matrix requires a scalar parameter");
        }
        return parameter;
    }

    private static int requireDimension(final int dimension) {
        if (dimension < 1) {
            throw new IllegalArgumentException("dimension must be at least 1");
        }
        return dimension;
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
            final Parameter scale = (Parameter) xo.getChild(Parameter.class);
            if (dimension < 1) {
                throw new XMLParseException(DIMENSION + " must be at least 1");
            }
            return new ScaledIdentityMatrixParameter(name, dimension, scale);
        }

        @Override
        public String getParserDescription() {
            return "A square matrix parameter equal to one scalar parameter times the identity.";
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return RULES;
        }

        @Override
        public Class getReturnType() {
            return MatrixParameter.class;
        }
    };

    private static final XMLSyntaxRule[] RULES = new XMLSyntaxRule[]{
            AttributeRule.newIntegerRule(DIMENSION, false),
            new ElementRule(Parameter.class, 1, 1)
    };
}
