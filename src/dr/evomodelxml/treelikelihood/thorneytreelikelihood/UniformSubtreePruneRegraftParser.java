package dr.evomodelxml.treelikelihood.thorneytreelikelihood;

import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.thorneytreelikelihood.ConstrainedTreeModel;
import dr.evomodel.treelikelihood.thorneytreelikelihood.ConstrainedTreeOperator;
import dr.evomodel.treelikelihood.thorneytreelikelihood.UniformSubtreePruneRegraft;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;


public class UniformSubtreePruneRegraftParser extends AbstractXMLObjectParser {
    public static final String UNIFORM_SUBTREE_PRUNE_REGRAFT = "uniformSubtreePruneRegraft";
    @Override
    public Object parseXMLObject(XMLObject  xo) throws XMLParseException {
        TreeModel tree =  (TreeModel) xo.getChild(TreeModel.class);
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);


       UniformSubtreePruneRegraft op = new UniformSubtreePruneRegraft(tree,weight);
       if(tree instanceof ConstrainedTreeModel){
           return new ConstrainedTreeOperator((ConstrainedTreeModel) tree, weight, op);
       }else{
           return op;
       }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return "An operator that prunes a regrafts within a clade.";
    }

    @Override
    public Class getReturnType() {
        return UniformSubtreePruneRegraft.class;
    }

    @Override
    public String getParserName() {
        return UNIFORM_SUBTREE_PRUNE_REGRAFT;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            new ElementRule(TreeModel.class)
    };
}
