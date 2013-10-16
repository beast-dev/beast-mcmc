package dr.evomodelxml.operators;

import dr.evomodel.operators.RateExchangeOperator;
import dr.evomodel.operators.SubtreeSlideOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class RateExchangeOperatorParser extends AbstractXMLObjectParser {

    public static final String RATE_EXCHANGE = "rateExchange";
    public static final String SWAP_TRAITS = "swapTraits";
    public static final String SWAP_RATES = "swapRates";
    public static final String SWAP_AT_ROOT = "swapAtRoot";
    public static final String MOVE_HEIGHT = "moveHeight";

    public String getParserName() {
        return RATE_EXCHANGE;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        boolean swapRates = xo.getBooleanAttribute(SWAP_RATES);
        boolean swapTraits = xo.getBooleanAttribute(SWAP_TRAITS);
        boolean swapAtRoot = xo.getBooleanAttribute(SWAP_AT_ROOT);
        boolean moveHeight = xo.getBooleanAttribute(MOVE_HEIGHT);
        return new RateExchangeOperator(treeModel, weight, swapRates, swapTraits, swapAtRoot, moveHeight);
    }

    public String getParserDescription() {
        return "An operator that exchanges rates and traits on a tree.";
    }

    public Class getReturnType() {
        return SubtreeSlideOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newBooleanRule(SWAP_RATES),
            AttributeRule.newBooleanRule(SWAP_TRAITS),
            AttributeRule.newBooleanRule(SWAP_AT_ROOT),
            AttributeRule.newBooleanRule(MOVE_HEIGHT),
            new ElementRule(TreeModel.class)
    };
}
