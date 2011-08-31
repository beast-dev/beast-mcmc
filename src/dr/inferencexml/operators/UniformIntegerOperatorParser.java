package dr.inferencexml.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.UniformIntegerOperator;
import dr.xml.*;

/**
 */
public class UniformIntegerOperatorParser extends AbstractXMLObjectParser {

    public final static String UNIFORM_INTEGER_OPERATOR = "uniformIntegerOperator";

    public String getParserName() {
        return UNIFORM_INTEGER_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);

        int lower = (int) (double) parameter.getBounds().getLowerLimit(0);
        if (xo.hasAttribute("lower")) lower = xo.getIntegerAttribute("lower");

        int upper = (int) (double) parameter.getBounds().getUpperLimit(0);
        if (xo.hasAttribute("upper")) upper = xo.getIntegerAttribute("upper");

        int count = 1;
        if (xo.hasAttribute("count")) count = xo.getIntegerAttribute("count");

        if (upper == lower || lower == (int) Double.NEGATIVE_INFINITY || upper == (int) Double.POSITIVE_INFINITY) {
            throw new XMLParseException(this.getParserName() + " boundaries not found in parameter "
                    + parameter.getParameterName() + " Use operator lower and upper !");
        }

        return new UniformIntegerOperator(parameter, lower, upper, weight, count);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "An operator that picks new parameter values uniformly at random.";
    }

    public Class getReturnType() {
        return UniformIntegerOperator.class;
    }


    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newDoubleRule("upper", true),
            AttributeRule.newDoubleRule("lower", true),
            AttributeRule.newDoubleRule("count", true),
            new ElementRule(Parameter.class)
    };
}
