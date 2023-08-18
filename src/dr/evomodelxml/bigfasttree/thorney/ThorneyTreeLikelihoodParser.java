
package dr.evomodelxml.bigfasttree.thorney;

import dr.evolution.datatype.ContinuousDataType;
import dr.evolution.datatype.IntegerDataType;
import dr.evomodel.bigfasttree.thorney.BranchLengthLikelihoodDelegate;
import dr.evomodel.bigfasttree.thorney.MutationBranchMap;
import dr.evomodel.bigfasttree.thorney.PoissonBranchLengthLikelihoodDelegate;
import dr.evomodel.bigfasttree.thorney.ThorneyDataLikelihoodDelegate;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.Likelihood;
import dr.xml.*;

/**
 * @author Andrew Rambaut
 * @author JT McCrone
 * @version $Id$
 */
public class ThorneyTreeLikelihoodParser extends AbstractXMLObjectParser {

    public static final String THORNEY_DATA_LIKELIHOOD_DELEGATE = "thorneyTreeLikelihood";

    public String getParserName() {
        return THORNEY_DATA_LIKELIHOOD_DELEGATE;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

        MutationBranchMap branchLengthProvider = (MutationBranchMap) xo.getChild(MutationBranchMap.class);
        BranchLengthLikelihoodDelegate branchLengthLikelihoodDelegate = (BranchLengthLikelihoodDelegate) xo.getChild(BranchLengthLikelihoodDelegate.class);
        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);
        ThorneyDataLikelihoodDelegate  dataLikelihoodDelegate = new ThorneyDataLikelihoodDelegate( treeModel, branchLengthProvider, branchLengthLikelihoodDelegate);
        
        return new TreeDataLikelihood(dataLikelihoodDelegate, treeModel, branchRateModel);
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
            new ElementRule(MutationBranchMap.class),
            new ElementRule(BranchLengthLikelihoodDelegate.class),
            new ElementRule(BranchRateModel.class)
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}


