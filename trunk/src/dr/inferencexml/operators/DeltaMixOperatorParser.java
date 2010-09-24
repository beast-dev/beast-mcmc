package dr.inferencexml.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.CoercableMCMCOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.DeltaMixOperator;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class DeltaMixOperatorParser extends AbstractXMLObjectParser {

    public static final String DELTA_MIX = "deltaMixOperator";
    public static final String DELTA = "delta";
    public static final String PARAMETER_WEIGHTS = "parameterWeights";

    public String getParserName() {
        return DELTA_MIX;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        CoercionMode mode = CoercionMode.parseMode(xo);

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        double delta = xo.getDoubleAttribute(DELTA);

        if (delta <= 0.0) {
            throw new XMLParseException("delta must be greater than 0.0");
        }

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);


        int[] parameterWeights;
        if (xo.hasAttribute(PARAMETER_WEIGHTS)) {
            parameterWeights = xo.getIntegerArrayAttribute(PARAMETER_WEIGHTS);
            System.out.print("Parameter weights for delta exchange are: ");
            for (int parameterWeight : parameterWeights) {
                System.out.print(parameterWeight + "\t");
            }
            System.out.println();

        } else {
            parameterWeights = new int[parameter.getDimension()];
            for (int i = 0; i < parameterWeights.length; i++) {
                parameterWeights[i] = 1;
            }
        }

        if (parameterWeights.length != parameter.getDimension()) {
            throw new XMLParseException("parameter weights have the same length as parameter");
        }


        return new DeltaMixOperator(parameter, parameterWeights, delta, weight, mode);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a scale operator on a given parameter.";
    }

    public Class getReturnType() {
        return MCMCOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(DELTA),
            AttributeRule.newIntegerArrayRule(PARAMETER_WEIGHTS, true),
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newBooleanRule(CoercableMCMCOperator.AUTO_OPTIMIZE, true),
            new ElementRule(Parameter.class)
    };
}