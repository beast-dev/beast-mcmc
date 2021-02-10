
package dr.evomodelxml.treelikelihood.thorneytreelikelihood;

import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.thorneytreelikelihood.BranchLengthProvider;
import dr.evomodel.treelikelihood.thorneytreelikelihood.ThorneyBranchLengthLikelihoodDelegate;
import dr.evomodel.treelikelihood.thorneytreelikelihood.ThorneyTreeLikelihood;
import dr.inference.model.Likelihood;
import dr.xml.*;

/**
 * @author Andrew Rambaut
 * @author JT McCrone
 * @version $Id$
 */
public class ThorneyTreeLikelihoodParser extends AbstractXMLObjectParser {

    public static final String TREE_LIKELIHOOD = "thorneyTreeLikelihood";

    public String getParserName() {
        return TREE_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class); // can be null
        BranchLengthProvider branchLengthProvider = (BranchLengthProvider) xo.getChild(BranchLengthProvider.class);
        ThorneyBranchLengthLikelihoodDelegate thorneyBranchLengthLikelihoodDelegate = (ThorneyBranchLengthLikelihoodDelegate) xo.getChild(ThorneyBranchLengthLikelihoodDelegate.class);
        return new ThorneyTreeLikelihood(TREE_LIKELIHOOD, treeModel,branchLengthProvider, thorneyBranchLengthLikelihoodDelegate, branchRateModel);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of a tree given  the number of mutations along each branch.";
    }

    public Class getReturnType() {
        return Likelihood.class;
    }

    public static final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class),
            new ElementRule(BranchLengthProvider.class),
            new ElementRule(ThorneyBranchLengthLikelihoodDelegate.class),
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}


