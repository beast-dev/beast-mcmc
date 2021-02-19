package dr.inferencexml.operators;

import dr.inference.operators.RepeatOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.xml.*;

public class RepeatOperatorParser extends AbstractXMLObjectParser {

    public static final String REPEAT_OPERATOR = "repeatOperator";
    public static final String WEIGHT = "weight";
    public static final String TARGET_ACCEPTANCE = "targetAcceptance";
    public static final String NONSHIFTED_MEAN = "nonshiftedMean";

    public String getParserName() {
        return REPEAT_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final double weight = xo.getDoubleAttribute(WEIGHT);

        final double targetProb = xo.getAttribute(TARGET_ACCEPTANCE, 0.234);
        final double nonshiftedMean = xo.getAttribute(NONSHIFTED_MEAN, 0.001);

        final SimpleMCMCOperator simpleOperator = (SimpleMCMCOperator) xo.getChild(0);

        if (targetProb <= 0.0 || targetProb >= 1.0)
            throw new RuntimeException("Target acceptance probability must be between 0.0 and 1.0");

        return new RepeatOperator(weight, targetProb,simpleOperator,nonshiftedMean);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element is a tunable operator that calls a nontunable operator an adaptable numbers of time per move";
    }

    public Class getReturnType() {
        return RepeatOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(SimpleMCMCOperator.class),
            AttributeRule.newDoubleRule(WEIGHT),
            AttributeRule.newDoubleRule(TARGET_ACCEPTANCE, true),
            AttributeRule.newDoubleRule(NONSHIFTED_MEAN, true)
    };
}