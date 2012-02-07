package dr.inferencexml.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.CoercableMCMCOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.ScaleOperator;
import dr.xml.*;

/**
 */
public class ScaleOperatorParser extends AbstractXMLObjectParser {
    public static final String SCALE_OPERATOR = "scaleOperator";
    public static final String SCALE_ALL = "scaleAll";
    public static final String SCALE_ALL_IND = "scaleAllIndependently";
    public static final String SCALE_FACTOR = "scaleFactor";
    public static final String DEGREES_OF_FREEDOM = "df";
    public static final String INDICATORS = "indicators";
    public static final String PICKONEPROB = "pickoneprob";

    public String getParserName() {
        return SCALE_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final boolean scaleAll = xo.getAttribute(SCALE_ALL, false);
        final boolean scaleAllInd = xo.getAttribute(SCALE_ALL_IND, false);
        final int degreesOfFreedom = xo.getAttribute(DEGREES_OF_FREEDOM, 0);

        final CoercionMode mode = CoercionMode.parseMode(xo);

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);

        if (scaleFactor <= 0.0 || scaleFactor >= 1.0) {
            throw new XMLParseException("scaleFactor must be between 0.0 and 1.0");
        }

        final Parameter parameter = (Parameter) xo.getChild(Parameter.class);

        Parameter indicator = null;
        double indicatorOnProb = 1.0;
        final XMLObject inds = xo.getChild(INDICATORS);

        if (inds != null) {
            indicator = (Parameter) inds.getChild(Parameter.class);
            if (inds.hasAttribute(PICKONEPROB)) {
                indicatorOnProb = inds.getDoubleAttribute(PICKONEPROB);
                if (!(0 <= indicatorOnProb && indicatorOnProb <= 1)) {
                    throw new XMLParseException("pickoneprob must be between 0.0 and 1.0");
                }
            }
        }
        ScaleOperator operator = new ScaleOperator(parameter, scaleAll,
                degreesOfFreedom, scaleFactor,
                mode, indicator, indicatorOnProb,
                scaleAllInd);
        operator.setWeight(weight);
        return operator;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a scale operator on a given parameter.";
    }

    public Class getReturnType() {
        return ScaleOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(SCALE_FACTOR),
            AttributeRule.newBooleanRule(SCALE_ALL, true),
            AttributeRule.newBooleanRule(SCALE_ALL_IND, true),
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newBooleanRule(CoercableMCMCOperator.AUTO_OPTIMIZE, true),
            AttributeRule.newIntegerRule(DEGREES_OF_FREEDOM, true),

            new ElementRule(Parameter.class),
            new ElementRule(INDICATORS,
                    new XMLSyntaxRule[]{
                            AttributeRule.newDoubleRule(PICKONEPROB, true),
                            new ElementRule(Parameter.class)}, true),
    };
}
