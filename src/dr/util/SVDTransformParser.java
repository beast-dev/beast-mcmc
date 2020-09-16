package dr.util;

import dr.xml.*;

public class SVDTransformParser extends AbstractXMLObjectParser {

    private static final String SVD_TRANSFORM = "svdTransform";
    private static final String NCOLS = "nColumns";
    private static final String NROWS = "nRows";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        int nRows = xo.getIntegerAttribute(NROWS);
        int nCols = xo.getIntegerAttribute(NCOLS);
        return new SVDTransform(nCols, nRows);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newIntegerRule(NCOLS),
                AttributeRule.newIntegerRule(NROWS)
        };
    }

    @Override
    public String getParserDescription() {
        return "Rotates a matrix to its orthogonal subspace.";
    }

    @Override
    public Class getReturnType() {
        return SVDTransform.class;
    }

    @Override
    public String getParserName() {
        return SVD_TRANSFORM;
    }
}
