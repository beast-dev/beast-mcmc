package dr.inferencexml.operators;

import dr.inference.model.GeneralParameterBounds;
import dr.inference.operators.SimpleMCMCOperator;
import dr.inference.operators.TransformedParameterOperator;
import dr.xml.*;

public class TransformedParameterOperatorParser extends AbstractXMLObjectParser {

    private static final String TRANSFORMED_OPERATOR = "transformedParameterOperator";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        SimpleMCMCOperator operator = (SimpleMCMCOperator) xo.getChild(SimpleMCMCOperator.class);
        GeneralParameterBounds bounds = (GeneralParameterBounds) xo.getChild(GeneralParameterBounds.class);
        return new TransformedParameterOperator(operator, bounds);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(SimpleMCMCOperator.class),
                new ElementRule(GeneralParameterBounds.class, true)
        };
    }

    @Override
    public String getParserDescription() {
        return "operator that corrects the hastings ratio with appropriate Jacobian term due to parameter transform";
    }

    @Override
    public Class getReturnType() {
        return TransformedParameterOperator.class;
    }

    @Override
    public String getParserName() {
        return TRANSFORMED_OPERATOR;
    }
}
