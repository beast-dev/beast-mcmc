package dr.evomodelxml.bigfasttree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.bigfasttree.ApproximatePoissonTreeLikelihood;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.xml.*;

import java.util.ArrayList;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class ApproximatePoissonTreeLikelihoodParser extends AbstractXMLObjectParser {

    public static final String TREE_LIKELIHOOD = "approximatePoissonTreeLikelihood";
    public static final String DATA = "data";

    public String getParserName() {
        return TREE_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        int sequenceLength = xo.getIntegerAttribute("sequenceLength");
        Tree dataTree = (Tree) xo.getElementFirstChild("data");
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

        checkTreesMatch(dataTree, treeModel);

        return new ApproximatePoissonTreeLikelihood(TREE_LIKELIHOOD, dataTree, sequenceLength, treeModel, branchRateModel);
    }

    /**
     * Helper function used in checkTreesMatch
     * @param tree
     * @return ArrayList of taxon ids
     */
    private static ArrayList<String> getTaxon(Tree tree) {
        ArrayList<String> taxa= new ArrayList<>();
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            NodeRef n = tree.getExternalNode(i);
            taxa.add((tree.getNodeTaxon(n).getId()));
        }
        return taxa;
    }

    /**
     * A helper function that ensures the data tree and tree model have the same taxa
     * @param dataTree
     * @param tree
     * @throws XMLParseException
     */
    private static void checkTreesMatch(Tree dataTree, Tree tree) throws XMLParseException {
        ArrayList<String> referenceTaxon = getTaxon(dataTree);
        ArrayList<String> treeTaxon = getTaxon(tree);

        if(treeTaxon.size()!=referenceTaxon.size()){
            throw new XMLParseException("TreeModel and data tree must have the same taxa");
        }

        for (String taxonName :
                treeTaxon) {
            if (!referenceTaxon.contains(taxonName)) {
                throw new XMLParseException("TreeModel and data tree must have the same taxa");
            }
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of a patternlist on a tree given the site model.";
    }

    public Class getReturnType() {
        return Likelihood.class;
    }

    public static final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule("sequenceLength", false),
            new ElementRule(DATA, new XMLSyntaxRule[] {
                    new ElementRule(Tree.class)
            }),
            new ElementRule(TreeModel.class),
            new ElementRule(BranchRateModel.class, true),
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}
