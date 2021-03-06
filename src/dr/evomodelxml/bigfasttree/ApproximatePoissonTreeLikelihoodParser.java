
package dr.evomodelxml.bigfasttree;

import dr.evolution.distance.DistanceMatrix;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.bigfasttree.*;
import dr.evomodel.bigfasttree.constrainedtree.CladeNodeModel;
import dr.evomodel.bigfasttree.constrainedtree.ConstrainedTreeBranchLengthProvider;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.xml.*;

/**
 * @author Andrew Rambaut
 * @author JT McCrone
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

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        BranchLengthProvider branchLengthProvider;
        if(xo.getElementFirstChild(DATA) instanceof  Tree){
            Tree dataTree = (Tree) xo.getElementFirstChild(DATA);
            CladeNodeModel cladeNodeModel = null;
            try {
                cladeNodeModel = new CladeNodeModel(dataTree, treeModel);
            } catch (TreeUtils.MissingTaxonException e) {
                e.printStackTrace();
            }
            branchLengthProvider = new ConstrainedTreeBranchLengthProvider(cladeNodeModel);
        }else if(xo.getElementFirstChild(DATA) instanceof  CladeNodeModel) {
            CladeNodeModel cladeNodeModel = (CladeNodeModel) xo.getElementFirstChild(DATA);
            branchLengthProvider = new ConstrainedTreeBranchLengthProvider(cladeNodeModel);
        }else{
                DistanceMatrix dataMatrix = (DistanceMatrix)xo.getElementFirstChild(DATA);
                branchLengthProvider = new RzhetskyNeiBranchLengthProvider(dataMatrix,treeModel);
            }

        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

        return new ApproximatePoissonTreeLikelihood(TREE_LIKELIHOOD, sequenceLength, treeModel, branchRateModel,branchLengthProvider);
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
//            new OrRule(
//                    new ElementRule(DATA, new XMLSyntaxRule[]{
//                            new ElementRule(Tree.class)
//                    }),
//                    new ElementRule(DATA, new XMLSyntaxRule[]{
//                            new ElementRule(DistanceMatrix.class)
//                    }),
//            ),
            new ElementRule(TreeModel.class),
            new ElementRule(BranchRateModel.class, true),
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}
