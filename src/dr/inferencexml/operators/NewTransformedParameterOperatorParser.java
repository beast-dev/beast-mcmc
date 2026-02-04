package dr.inferencexml.operators;

import dr.inference.model.BoundedSpace;
import dr.inference.operators.NewTransformedParameterOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.util.Transform;
import dr.xml.*;

import static dr.inference.operators.NewTransformedParameterOperator.NEW_TRANSFORMED_OPERATOR;

public class NewTransformedParameterOperatorParser extends AbstractXMLObjectParser {

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        SimpleMCMCOperator operator = (SimpleMCMCOperator) xo.getChild(SimpleMCMCOperator.class);
        BoundedSpace bounds = (BoundedSpace) xo.getChild(BoundedSpace.class);
        Transform transform = (Transform) xo.getChild(Transform.class);
        return new NewTransformedParameterOperator(operator, transform, bounds);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(SimpleMCMCOperator.class),
                new ElementRule(Transform.class),
                new ElementRule(BoundedSpace.class, true)
        };
    }

    @Override
    public String getParserDescription() {
        return "Operator that corrects the hastings ratio with appropriate Jacobian term due to parameter transform";
    }

    @Override
    public Class getReturnType() {
        return NewTransformedParameterOperator.class;
    }

    @Override
    public String getParserName() {
        return NEW_TRANSFORMED_OPERATOR;
    }
}
