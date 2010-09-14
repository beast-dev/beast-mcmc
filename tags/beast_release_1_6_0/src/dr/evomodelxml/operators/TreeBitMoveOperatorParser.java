package dr.evomodelxml.operators;

import dr.evomodel.operators.TreeBitMoveOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class TreeBitMoveOperatorParser extends AbstractXMLObjectParser {

    public static final String BIT_MOVE_OPERATOR = "treeBitMoveOperator";
    public static final String INDICTATOR_TRAIT = "indicatorTrait";
    public static final String TRAIT2 = "trait2";

    public String getParserName() {
        return BIT_MOVE_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);


        String trait1 = null;
        String trait2 = null;
        if (xo.hasAttribute(INDICTATOR_TRAIT)) trait1 = xo.getStringAttribute(INDICTATOR_TRAIT);
        if (xo.hasAttribute(TRAIT2)) trait2 = xo.getStringAttribute(TRAIT2);

        return new TreeBitMoveOperator(treeModel, trait1, trait2, weight);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a bit-move operator on a given parameter.";
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
            AttributeRule.newStringRule(TRAIT2, true)
    };

}
