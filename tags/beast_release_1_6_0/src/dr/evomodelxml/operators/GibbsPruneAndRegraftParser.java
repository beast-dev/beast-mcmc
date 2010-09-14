package dr.evomodelxml.operators;

import dr.evomodel.operators.GibbsPruneAndRegraft;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class GibbsPruneAndRegraftParser extends AbstractXMLObjectParser {

    public static final String GIBBS_PRUNE_AND_REGRAFT = "GibbsPruneAndRegraft";

    public String getParserName() {
        return GIBBS_PRUNE_AND_REGRAFT;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        boolean pruned = true;
        if (xo.hasAttribute("pruned")) {
            pruned = xo.getBooleanAttribute("pruned");
        }

        return new GibbsPruneAndRegraft(treeModel, pruned, weight);
    }

    // ************************************************************************
    // AbstractXMLObjectParser implementation
    // ************************************************************************

    public String getParserDescription() {
        return "This element represents a Gibbs sampler implemented through a prune and regraft operator. "
                + "This operator prunes a random subtree and regrafts it below a node chosen by an importance distribution which is the proportion of the likelihoods of the proposals.";
    }

    public Class getReturnType() {
        return GibbsPruneAndRegraft.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules;

    {
        rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                AttributeRule.newBooleanRule("pruned"),
                new ElementRule(TreeModel.class)};
    }

}
