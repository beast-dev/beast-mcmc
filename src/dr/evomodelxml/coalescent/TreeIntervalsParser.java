package dr.evomodelxml.coalescent;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.TaxonList;
import dr.evomodel.coalescent.TreeIntervals;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

public class TreeIntervalsParser extends AbstractXMLObjectParser{

    public static final String TREE_INTERVALS = "treeIntervals";
    public static final String INCLUDE = "include";
    public static final String EXCLUDE = "exclude";

    public String getParserName() {
        return TREE_INTERVALS;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Tree tree = (Tree) xo.getChild(Tree.class);

        TaxonList includeSubtree = null;

        if (xo.hasChildNamed(INCLUDE)) {
            includeSubtree = (TaxonList) xo.getElementFirstChild(INCLUDE);
        }

        List<TaxonList> excludeSubtrees = new ArrayList<>();

        if (xo.hasChildNamed(EXCLUDE)) {
            XMLObject cxo = xo.getChild(EXCLUDE);
            for (int i = 0; i < cxo.getChildCount(); i++) {
                excludeSubtrees.add((TaxonList) cxo.getChild(i));
            }
        }

        try {
            return new TreeIntervals(tree, includeSubtree, excludeSubtrees);
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