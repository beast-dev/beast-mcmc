package dr.inference.model;

import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Backward-compatible XML alias for column-wise {@link NormalizedMatrixParameter}.
 */
public final class ColumnNormalizedMatrixParameter extends NormalizedMatrixParameter {

    public static final String NAME = "columnNormalizedMatrixParameter";

    private static final String MIN_COLUMN_NORM = "minColumnNorm";

    public ColumnNormalizedMatrixParameter(final String name,
                                           final MatrixParameter rawMatrix) {
        this(name, rawMatrix, DEFAULT_MIN_NORM);
    }

    public ColumnNormalizedMatrixParameter(final String name,
                                           final MatrixParameter rawMatrix,
                                           final double minColumnNorm) {
        super(name, rawMatrix, Normalization.COLUMNS, minColumnNorm);
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
            final double minColumnNorm = xo.getAttribute(MIN_COLUMN_NORM, DEFAULT_MIN_NORM);
            return new ColumnNormalizedMatrixParameter(name, raw, minColumnNorm);
        }

        @Override
        public String getParserDescription() {
            return "A dense matrix parameter whose columns are normalized versions of a raw dense matrix.";
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return RULES;
        }

        @Override
        public Class getReturnType() {
            return ColumnNormalizedMatrixParameter.class;
        }
    };

    private static final XMLSyntaxRule[] RULES = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(MIN_COLUMN_NORM, true),
            new ElementRule(MatrixParameter.class)
    };
}
