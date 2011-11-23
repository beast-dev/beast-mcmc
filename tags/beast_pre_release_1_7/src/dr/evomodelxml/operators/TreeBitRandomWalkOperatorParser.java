package dr.evomodelxml.operators;

import dr.evomodel.operators.TreeBitRandomWalkOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class TreeBitRandomWalkOperatorParser extends AbstractXMLObjectParser {

    public static final String BIT_RANDOM_WALK_OPERATOR = "treeBitRandomWalk";
    public static final String INDICTATOR_TRAIT = "indicatorTrait";
    public static final String TRAIT2 = "trait2";
    public static final String SWAP_TRAIT2 = "swapTrait2";

    public String getParserName() {
        return BIT_RANDOM_WALK_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);


        String trait1 = null;
        String trait2 = null;
        if (xo.hasAttribute(INDICTATOR_TRAIT)) trait1 = xo.getStringAttribute(INDICTATOR_TRAIT);
        if (xo.hasAttribute(TRAIT2)) trait2 = xo.getStringAttribute(TRAIT2);
        int k = xo.getAttribute("k", 1);
        boolean swapTrait2 = xo.getAttribute(SWAP_TRAIT2, true);

        return new TreeBitRandomWalkOperator(treeModel, trait1, trait2, weight, k, swapTrait2);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a bit-random walk operator on a random " +
                "indicator/variable pair in the tree.";
    }

    public Class getReturnType() {
        return MCMCOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            new ElementRule(TreeModel.class),
            AttributeRule.newStringRule(INDICTATOR_TRAIT, true),
            AttributeRule.newStringRule(TRAIT2, true),
            AttributeRule.newBooleanRule(SWAP_TRAIT2, true),
            AttributeRule.newIntegerRule("k", true)
    };
}
