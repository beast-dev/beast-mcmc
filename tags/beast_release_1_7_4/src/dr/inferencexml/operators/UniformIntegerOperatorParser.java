package dr.inferencexml.operators;

import dr.inference.model.Parameter;
import dr.inference.model.Variable;
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

        Variable parameter = (Variable) xo.getChild(Variable.class);

        int count = 1;
        if (xo.hasAttribute("count")) count = xo.getIntegerAttribute("count");

        if (parameter instanceof Parameter) {
            int lower = (int) (double) ((Parameter) parameter).getBounds().getLowerLimit(0);
            if (xo.hasAttribute("lower")) lower = xo.getIntegerAttribute("lower");

            int upper = (int) (double) ((Parameter) parameter).getBounds().getUpperLimit(0);
            if (xo.hasAttribute("upper")) upper = xo.getIntegerAttribute("upper");

            if (upper == lower || lower == (int) Double.NEGATIVE_INFINITY || upper == (int) Double.POSITIVE_INFINITY) {
                throw new XMLParseException(this.getParserName() + " boundaries not found in parameter "
                        + parameter.getId() + " Use operator lower and upper !");
            }

            return new UniformIntegerOperator((Parameter) parameter, lower, upper, weight, count);
        } else { // Variable<Integer>, Bounds.Staircase
            return new UniformIntegerOperator(parameter, weight, count);
        }
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
            new ElementRule(Variable.class)
    };
}
