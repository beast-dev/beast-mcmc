package dr.evomodelxml.operators;

import dr.evomodel.operators.SubtreeJumpOperator;
import dr.evomodel.operators.SubtreeSlideOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.CoercableMCMCOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class SubtreeJumpOperatorParser extends AbstractXMLObjectParser {

    public static final String SUBTREE_JUMP = "subtreeJump";

    public String getParserName() {
        return SUBTREE_JUMP;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        /*
        No coercion at the moment.
        CoercionMode mode = CoercionMode.DEFAULT;
        if (xo.hasAttribute(CoercableMCMCOperator.AUTO_OPTIMIZE)) {
            if (xo.getBooleanAttribute(CoercableMCMCOperator.AUTO_OPTIMIZE)) {
                mode = CoercionMode.COERCION_ON;
            } else {
                mode = CoercionMode.COERCION_OFF;
            }
        }
        */

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

//        final double targetAcceptance = xo.getAttribute(TARGET_ACCEPTANCE, 0.234);
//
//        final double size = xo.getAttribute("size", 1.0);
//
//        if (Double.isInfinite(size) || size <= 0.0) {
//            throw new XMLParseException("size attribute must be positive and not infinite. was " + size +
//           " for tree " + treeModel.getId() );
//        }

        SubtreeJumpOperator operator = new SubtreeJumpOperator(treeModel, weight);
//        operator.setTargetAcceptanceProbability(targetAcceptance);

        return operator;
    }

    public String getParserDescription() {
        return "An operator that jumps a subtree to another edge at the same height.";
    }

    public Class getReturnType() {
        return SubtreeJumpOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            // No coercion at the moment.
            //AttributeRule.newDoubleRule("size", true),
            //AttributeRule.newBooleanRule(CoercableMCMCOperator.AUTO_OPTIMIZE, true),
            new ElementRule(TreeModel.class)
    };

}
