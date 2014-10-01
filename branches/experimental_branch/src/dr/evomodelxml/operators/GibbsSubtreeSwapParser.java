package dr.evomodelxml.operators;

import dr.evomodel.operators.GibbsSubtreeSwap;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class GibbsSubtreeSwapParser extends AbstractXMLObjectParser {

    public static final String GIBBS_SUBTREE_SWAP = "GibbsSubtreeSwap";

    public String getParserName() {
        return GIBBS_SUBTREE_SWAP;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final boolean pruned = xo.getAttribute("pruned", true);

        return new GibbsSubtreeSwap(treeModel, pruned, weight);
    }

    // ************************************************************************
    // AbstractXMLObjectParser implementation
    // ************************************************************************

    public String getParserDescription() {
        return "This element represents a Gibbs wide exchange operator. "
                + "This operator swaps two subtrees chosen to their posterior probaility.";
    }

    public Class getReturnType() {
        return GibbsSubtreeSwap.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newBooleanRule("pruned", true),
            new ElementRule(TreeModel.class)
    };

}
