package dr.inferencexml.operators;

import dr.inference.model.BoundedSpace;
import dr.inference.operators.SimpleMCMCOperator;
import dr.inference.operators.TransformedParameterOperator;
import dr.xml.*;

import static dr.inference.operators.TransformedParameterOperator.TRANSFORMED_OPERATOR;

public class TransformedParameterOperatorParser extends AbstractXMLObjectParser {


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        SimpleMCMCOperator operator = (SimpleMCMCOperator) xo.getChild(SimpleMCMCOperator.class);
        BoundedSpace bounds = (BoundedSpace) xo.getChild(BoundedSpace.class);
        return new TransformedParameterOperator(operator, bounds);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(SimpleMCMCOperator.class),
                new ElementRule(BoundedSpace.class, true)
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
