package dr.inferencexml.operators;

import dr.inference.operators.JointOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.xml.*;

/**
 *
 */
public class JointOperatorParser extends AbstractXMLObjectParser {

    public static final String JOINT_OPERATOR = "jointOperator";
    public static final String WEIGHT = "weight";
    public static final String TARGET_ACCEPTANCE = "targetAcceptance";

    public String getParserName() {
        return JOINT_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final double weight = xo.getDoubleAttribute(WEIGHT);

        final double targetProb = xo.getAttribute(TARGET_ACCEPTANCE, 0.2);

        if (targetProb <= 0.0 || targetProb >= 1.0)
            throw new RuntimeException("Target acceptance probability must be between 0.0 and 1.0");

        JointOperator operator = new JointOperator(weight, targetProb);

        for (int i = 0; i < xo.getChildCount(); i++) {
            operator.addOperator((SimpleMCMCOperator) xo.getChild(i));
        }

        return operator;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents an arbitrary list of operators; only the first is optimizable";
    }

    public Class getReturnType() {
        return JointOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(SimpleMCMCOperator.class, 1, Integer.MAX_VALUE),
            AttributeRule.newDoubleRule(WEIGHT),
            AttributeRule.newDoubleRule(TARGET_ACCEPTANCE, true)
    };
}
