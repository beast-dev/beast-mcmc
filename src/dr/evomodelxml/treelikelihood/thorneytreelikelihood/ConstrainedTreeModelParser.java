package dr.evomodelxml.treelikelihood.thorneytreelikelihood;

import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.thorneytreelikelihood.ConstrainedTreeModel;

import dr.xml.*;

import java.util.logging.Logger;

public class ConstrainedTreeModelParser extends AbstractXMLObjectParser {

    public static final String CONSTRAINED_TREE_MODEL = "constrainedTreeModel";
    public static final String CONSTRAINTS_TREE = "constraintsTree";

    public String getParserName() {
        return CONSTRAINED_TREE_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Tree tree = (Tree) xo.getChild(Tree.class);
        Tree constraintsTree = (Tree) xo.getElementFirstChild(CONSTRAINTS_TREE);
        Logger.getLogger("dr.evomodel").info("\nCreating the constrained tree model based on big fast tree model, '" + xo.getId() + "'");
        ConstrainedTreeModel treeModel =  new ConstrainedTreeModel(xo.getId(),tree, constraintsTree);
        Logger.getLogger("dr.evomodel").info("  taxon count = " + treeModel.getExternalNodeCount());
        Logger.getLogger("dr.evomodel").info("  tree height = " + treeModel.getNodeHeight(treeModel.getRoot()));
        return treeModel;
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
            new ElementRule(Tree.class),
            new ElementRule(CONSTRAINTS_TREE, new XMLSyntaxRule[]{
                    new ElementRule(Tree.class)
            }),
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}
