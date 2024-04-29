package dr.evomodelxml.bigfasttree.thorney;

import dr.evomodel.bigfasttree.thorney.ConstrainedTreeModel;
import dr.evomodel.bigfasttree.thorney.ConstrainedTreeOperator;
import dr.evomodel.bigfasttree.thorney.UniformSubtreePruneRegraft;
import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.tree.TreeModel;
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
           return ConstrainedTreeOperator.parse((ConstrainedTreeModel) tree, weight, op,xo);
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
        return AbstractTreeOperator.class;
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
