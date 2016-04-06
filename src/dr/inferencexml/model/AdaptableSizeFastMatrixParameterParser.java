package dr.inferencexml.model;


import dr.xml.*;

/**
 * Created by max on 4/6/16.
 */
public class AdaptableSizeFastMatrixParameterParser extends AbstractXMLObjectParser {


    private static final String ADAPTABLE_SIZE_FAST_MATRIX_PARAMETER="adaptableSizeFastMatrixParameter";
    private static final String ROWS="rows";
    private static final String MAX_ROW_SIZE="maxRowSize";
    private static final String MAX_COL_SIZE="maxColumnSize";
    private static final String COLUMNS="columns";
    private static final String TRANSPOSE="transpose";
    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
      return 0;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
//            new ElementRule(Parameter.class, 0, Integer.MAX_VALUE),
            AttributeRule.newBooleanRule(TRANSPOSE, true),
            AttributeRule.newIntegerRule(COLUMNS, true),
            AttributeRule.newIntegerRule(ROWS, true),

    };


    @Override
    public String getParserDescription() {
        return "Returns a blockUpperTriangularMatrixParameter which is a compoundParameter which forces the last element to be of full length, the second to last element to be of full length-1, etc.";  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public Class getReturnType() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getParserName() {
        return ADAPTABLE_SIZE_FAST_MATRIX_PARAMETER;  //To change body of implemented methods use File | Settings | File Templates.
    }
}