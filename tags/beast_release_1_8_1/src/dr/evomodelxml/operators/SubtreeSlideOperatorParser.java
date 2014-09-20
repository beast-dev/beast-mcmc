package dr.evomodelxml.operators;

import dr.evomodel.operators.SubtreeSlideOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.CoercableMCMCOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class SubtreeSlideOperatorParser extends AbstractXMLObjectParser {

    public static final String SUBTREE_SLIDE = "subtreeSlide";
    public static final String SWAP_RATES = "swapInRandomRate";
    public static final String SWAP_TRAITS = "swapInRandomTrait";
    public static final String DIRICHLET_BRANCHES = "branchesAreScaledDirichlet";
    public static final String TARGET_ACCEPTANCE = "targetAcceptance";

    public static final String TRAIT = "trait";

    public String getParserName() {
        return SUBTREE_SLIDE;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean swapRates = xo.getAttribute(SWAP_RATES, false);
        boolean swapTraits = xo.getAttribute(SWAP_TRAITS, false);
        boolean scaledDirichletBranches = xo.getAttribute(DIRICHLET_BRANCHES, false);

        CoercionMode mode = CoercionMode.DEFAULT;
        if (xo.hasAttribute(CoercableMCMCOperator.AUTO_OPTIMIZE)) {
            if (xo.getBooleanAttribute(CoercableMCMCOperator.AUTO_OPTIMIZE)) {
                mode = CoercionMode.COERCION_ON;
            } else {
                mode = CoercionMode.COERCION_OFF;
            }
        }

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        final double targetAcceptance = xo.getAttribute(TARGET_ACCEPTANCE, 0.234);

        final double size = xo.getAttribute("size", 1.0);

        if (Double.isInfinite(size) || size <= 0.0) {
            throw new XMLParseException("size attribute must be positive and not infinite. was " + size +
           " for tree " + treeModel.getId() );
        }

        final boolean gaussian = xo.getBooleanAttribute("gaussian");
        SubtreeSlideOperator operator = new SubtreeSlideOperator(treeModel, weight, size, gaussian,
                swapRates, swapTraits, scaledDirichletBranches, mode);
        operator.setTargetAcceptanceProbability(targetAcceptance);

        return operator;
    }

    public String getParserDescription() {
        return "An operator that slides a subtree.";
    }

    public Class getReturnType() {
        return SubtreeSlideOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            // Make size optional. If not given or equals zero, size is set to half of average tree branch length.
            AttributeRule.newDoubleRule("size", true),
            AttributeRule.newDoubleRule(TARGET_ACCEPTANCE, true),
            AttributeRule.newBooleanRule("gaussian"),
            AttributeRule.newBooleanRule(SWAP_RATES, true),
            AttributeRule.newBooleanRule(SWAP_TRAITS, true),
            AttributeRule.newBooleanRule(CoercableMCMCOperator.AUTO_OPTIMIZE, true),
            new ElementRule(TreeModel.class)
    };

}
