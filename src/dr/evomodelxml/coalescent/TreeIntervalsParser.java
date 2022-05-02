package dr.evomodelxml.coalescent;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.coalescent.TreeIntervals;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

public class TreeIntervalsParser extends AbstractXMLObjectParser{

    public static final String TREE_INTERVALS = "treeIntervals";
    public static final String TREE = "tree";

    public String getParserName() {
        return TREE_INTERVALS;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Tree tree = (Tree) xo.getChild(Tree.class);

        try {
            return new TreeIntervals(tree, null, null);
        } catch (TreeUtils.MissingTaxonException mte) {
            throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + " was not found in the tree.");
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Parser for TreeIntervals.";
    }

    public Class getReturnType() {
        return TreeIntervals.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class)
    };

}