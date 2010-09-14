package dr.inferencexml.operators;

import dr.inference.operators.MCMCOperator;
import dr.inference.operators.TeamOperator;
import dr.xml.*;

/**
 *
 */
public class TeamOperatorParser extends AbstractXMLObjectParser {

    public static final String TEAM_OPERATOR = "teamOperator";
    public static final String SUBSET_SIZE = "size";
    // public static final String TARGET_ACCEPTANCE = "targetAcceptance";

    public String getParserName() {
        return TEAM_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final int size = xo.getIntegerAttribute(SUBSET_SIZE);
        //final double targetProb = xo.getAttribute(TARGET_ACCEPTANCE, 0.2);

//            if (targetProb <= 0.0 || targetProb >= 1.0)
//                throw new RuntimeException("Target acceptance probability must be between 0.0 and 1.0");

        MCMCOperator[] o = new MCMCOperator[xo.getChildCount()];
        for (int i = 0; i < o.length; i++) {
            o[i] = (MCMCOperator) xo.getChild(i);
        }

        return new TeamOperator(o, size, weight);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "An arbitrary list of operators; A random subset of size N is aggregated in one operation." +
                " Operators may have unequal weights - in that case a subset probability of selection is proportional to " +
                "the sum of it's members weights.";
    }

    public Class getReturnType() {
        return TeamOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(MCMCOperator.class, 1, Integer.MAX_VALUE),
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newIntegerRule(SUBSET_SIZE)
            //AttributeRule.newDoubleRule(TARGET_ACCEPTANCE, true)
    };
}
