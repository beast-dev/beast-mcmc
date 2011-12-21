package dr.inferencexml.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SetOperator;
import dr.xml.*;

/**
 */
public class SetOperatorParser extends AbstractXMLObjectParser {

    public static final String SET_OPERATOR = "setOperator";
    public static final String SET = "set";

    public String getParserName() {
        return SET_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double[] values = xo.getDoubleArrayAttribute(SET);
        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);

        System.out.println("Creating set operator for parameter " + parameter.getParameterName());
        System.out.print("  set = {" + values[0]);
        for (int i = 1; i < values.length; i++) {
            System.out.print(", " + values[i]);
        }
        System.out.println("}");

        SetOperator operator = new SetOperator(parameter, values);
        operator.setWeight(weight);

        return operator;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents an operator on a set.";
    }

    public Class getReturnType() {
        return SetOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleArrayRule(SET),
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            new ElementRule(Parameter.class)
    };
}
