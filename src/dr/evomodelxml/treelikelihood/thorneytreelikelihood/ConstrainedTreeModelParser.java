package dr.evomodelxml.treelikelihood.thorneytreelikelihood;

import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.thorneytreelikelihood.ConstrainedTreeModel;

import dr.xml.*;

public class ConstrainedTreeModelParser extends AbstractXMLObjectParser {

    public static final String CONSTRAINED_TREE_MODEL = "constrainedTreeModel";
    public static final String CONSTRAINTS_TREE = "constraintsTree";

    public String getParserName() {
        return CONSTRAINED_TREE_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        Tree constraintsTree = (Tree) xo.getElementFirstChild(CONSTRAINTS_TREE);
        return new ConstrainedTreeModel(treeModel, constraintsTree);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the a tree model that with defined clades that may not be broken.";
    }

    public Class getReturnType() {
        return TreeModel.class;
    }

    public static final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class),
            new ElementRule(CONSTRAINTS_TREE, new XMLSyntaxRule[]{
                    new ElementRule(Tree.class)
            }),
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}
