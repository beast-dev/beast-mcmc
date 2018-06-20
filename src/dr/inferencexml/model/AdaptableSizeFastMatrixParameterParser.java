package dr.inferencexml.model;


import dr.inference.model.AdaptableSizeFastMatrixParameter;
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
    private static final String STARTING_VALUE = "startingValue";
    private static final String LOWER_TRIANGLE = "lowerTriangle";
    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        int MaxRowSize = xo.getAttribute(MAX_ROW_SIZE, 1);
        int MaxColumnSize = xo.getAttribute(MAX_COL_SIZE, 1);
        int rowDimension = xo.getAttribute(ROWS, 1);
        int columnDimension = xo.getAttribute(COLUMNS, 1);
        String name = xo.getId();
        double startingValue = 1;
        if(xo.hasAttribute(STARTING_VALUE))
            startingValue = xo.getDoubleAttribute(STARTING_VALUE);
        boolean lowerTriangle = false;
        if(xo.hasAttribute(LOWER_TRIANGLE))
            lowerTriangle = xo.getBooleanAttribute(LOWER_TRIANGLE);

      return new AdaptableSizeFastMatrixParameter(name, rowDimension, columnDimension, MaxRowSize, MaxColumnSize, startingValue, lowerTriangle);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
//            new ElementRule(Parameter.class, 0, Integer.MAX_VALUE),
            AttributeRule.newIntegerRule(COLUMNS, true),
            AttributeRule.newIntegerRule(ROWS, true),
            AttributeRule.newIntegerRule(MAX_ROW_SIZE, true),
            AttributeRule.newIntegerRule(MAX_COL_SIZE, true),
            AttributeRule.newDoubleRule(STARTING_VALUE, true),
            AttributeRule.newBooleanRule(LOWER_TRIANGLE, true),

    };



    @Override
    public String getParserDescription() {
        return "Returns a blockUpperTriangularMatrixParameter which is a compoundParameter which forces the last element to be of full length, the second to last element to be of full length-1, etc.";  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public Class getReturnType() {
        return AdaptableSizeFastMatrixParameter.class;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getParserName() {
        return ADAPTABLE_SIZE_FAST_MATRIX_PARAMETER;  //To change body of implemented methods use File | Settings | File Templates.
    }
}