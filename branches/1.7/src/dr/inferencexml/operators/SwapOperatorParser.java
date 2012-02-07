package dr.inferencexml.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.SwapOperator;
import dr.xml.*;

/**
 */
public class SwapOperatorParser extends AbstractXMLObjectParser {

    public final static String SWAP_OPERATOR = "swapOperator";

    public String getParserName() {
        return SWAP_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        double weight = xo.getDoubleAttribute("weight");
        int size = xo.getIntegerAttribute("size");

        boolean autoOptimize = xo.getBooleanAttribute("autoOptimize");
        if (autoOptimize) throw new XMLParseException("swapOperator can't be optimized!");

        System.out.println("Creating swap operator for parameter " + parameter.getParameterName() + " (weight=" + weight + ")");

        SwapOperator so = new SwapOperator(parameter, size);
        so.setWeight(weight);

        return so;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents an operator that swaps values in a multi-dimensional parameter.";
    }

    public Class getReturnType() {
        return SwapOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{

            AttributeRule.newDoubleRule("weight"),
            AttributeRule.newIntegerRule("size"),
            AttributeRule.newBooleanRule("autoOptimize"),
            new ElementRule(Parameter.class)
    };

}
