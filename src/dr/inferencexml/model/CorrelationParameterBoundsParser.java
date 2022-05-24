package dr.inferencexml.model;

import dr.inference.model.GeneralParameterBounds;
import dr.xml.*;

public class CorrelationParameterBoundsParser extends AbstractXMLObjectParser {
    private static final String CORRELATION_BOUNDS = "correlationBounds";
    private static final String DIMENSION = "dimension";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        int dim = xo.getIntegerAttribute(DIMENSION);
        return new GeneralParameterBounds.CorrelationParameterBounds(dim);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newIntegerRule(DIMENSION)
        };
    }

    @Override
    public String getParserDescription() {
        return "Indicates whether a parameter is a valid correlation matrix or not";
    }

    @Override
    public Class getReturnType() {
        return GeneralParameterBounds.CorrelationParameterBounds.class;
    }

    @Override
    public String getParserName() {
        return CORRELATION_BOUNDS;
    }
}
