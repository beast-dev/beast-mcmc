package dr.evomodelxml.speciation;

import dr.evomodel.speciation.BranchingLikelihood;
import dr.evomodel.speciation.BranchingModel;
import dr.evomodel.tree.TreeModel;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLSyntaxRule;

/**
 */
public class BranchingLikelihoodParser extends AbstractXMLObjectParser {

    public static final String BRANCHING_LIKELIHOOD = "branchingLikelihood";
    public static final String MODEL = "model";
    public static final String TREE = "branchingTree";

    public String getParserName() {
        return BRANCHING_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) {

        XMLObject cxo = xo.getChild(MODEL);
        BranchingModel branchingModel = (BranchingModel) cxo.getChild(BranchingModel.class);

        cxo = xo.getChild(TREE);
        TreeModel treeModel = (TreeModel) cxo.getChild(TreeModel.class);

        return new BranchingLikelihood(treeModel, branchingModel);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of the tree given the demographic function.";
    }

    public Class getReturnType() {
        return BranchingLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(MODEL, new XMLSyntaxRule[]{
                    new ElementRule(BranchingModel.class)
            }),
            new ElementRule(TREE, new XMLSyntaxRule[]{
                    new ElementRule(TreeModel.class)
            }),
    };

}
