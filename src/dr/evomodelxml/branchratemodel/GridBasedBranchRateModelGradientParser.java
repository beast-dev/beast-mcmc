package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.GridBasedBranchRateModel;
import dr.evomodel.branchratemodel.GridBasedBranchRateModelGradient;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;

public class GridBasedBranchRateModelGradientParser extends AbstractXMLObjectParser {
    public final String PARSER_NAME = "gridBasedBranchRateModelGradient";
    private final String TRAIT_NAME = "traitName";
    private final String USE_HESSIAN = "useHessian";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String traitName = xo.getAttribute(TRAIT_NAME, DEFAULT_TRAIT_NAME);
        boolean useHessian = xo.getAttribute(USE_HESSIAN, false);
        Parameter rates = (Parameter) xo.getChild(Parameter.class);
        GridBasedBranchRateModel branchRateModel = (GridBasedBranchRateModel) xo.getChild(GridBasedBranchRateModel.class);
        TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
        BeagleDataLikelihoodDelegate likelihoodDelegate = (BeagleDataLikelihoodDelegate) treeDataLikelihood.getDataLikelihoodDelegate();
        return new GridBasedBranchRateModelGradient(traitName, treeDataLikelihood, likelihoodDelegate, rates, branchRateModel, useHessian);
    }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(TRAIT_NAME, true),
            AttributeRule.newBooleanRule(USE_HESSIAN, true),
            new ElementRule(GridBasedBranchRateModel.class),
            new ElementRule(TreeDataLikelihood.class),
            new ElementRule(Parameter.class)
    };

    @Override
    public String getParserDescription() {
        return "Parser for grid based branch rate model";
    }

    @Override
    public Class getReturnType() {
        return GridBasedBranchRateModelGradient.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }
}
