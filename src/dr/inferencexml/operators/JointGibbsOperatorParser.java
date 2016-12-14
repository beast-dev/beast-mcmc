package dr.inferencexml.operators;

import dr.inference.operators.GibbsOperator;
import dr.inference.operators.JointGibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.xml.*;

/**
 * @author Max Tolkoff
 */
public class JointGibbsOperatorParser extends AbstractXMLObjectParser {

    public static final String JOINT_GIBBS_OPERATOR = "jointGibbsOperator";
    public static final String WEIGHT = "weight";
    public static final String TARGET_ACCEPTANCE = "targetAcceptance";

    public String getParserName() {
        return JOINT_GIBBS_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final double weight = xo.getDoubleAttribute(WEIGHT);

        JointGibbsOperator operator = new JointGibbsOperator(weight);

        for (int i = 0; i < xo.getChildCount(); i++) {
            if (xo.getChild(i) instanceof  GibbsOperator)
                operator.addOperator((SimpleMCMCOperator) xo.getChild(i));
            else
                throw new RuntimeException("Operator list must consist only of GibbsOperators");
        }

        return operator;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents a list of Gibbs Operators";
    }

    public Class getReturnType() {
        return JointGibbsOperator.class;
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
