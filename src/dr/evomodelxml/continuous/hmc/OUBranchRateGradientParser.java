package dr.evomodelxml.continuous.hmc;

import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DifferentiableBranchRates;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.OUBranchRateGradient;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.Parameter;
import dr.xml.*;

import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;

public final class OUBranchRateGradientParser extends AbstractXMLObjectParser {

    public static final String NAME = "ouBranchRateGradient";
    public static final String TRAIT_NAME = TreeTraitParserUtilities.TRAIT_NAME;

    @Override
    public String getParserName() {
        return NAME;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final String traitName = xo.getAttribute(TRAIT_NAME, DEFAULT_TRAIT_NAME);
        final TreeDataLikelihood treeDataLikelihood =
                (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

        final BranchRateModel branchRateModel = treeDataLikelihood.getBranchRateModel();
        if (!(branchRateModel instanceof DifferentiableBranchRates)) {
            throw new XMLParseException(NAME + " requires differentiable branch rates.");
        }

        final DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        if (!(delegate instanceof ContinuousDataLikelihoodDelegate)) {
            throw new XMLParseException(NAME + " requires a continuous traitDataLikelihood.");
        }

        final Parameter branchRates = ((DifferentiableBranchRates) branchRateModel).getRateParameter();
        return new OUBranchRateGradient(
                traitName,
                treeDataLikelihood,
                (ContinuousDataLikelihoodDelegate) delegate,
                branchRates);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(TRAIT_NAME, true),
            new ElementRule(TreeDataLikelihood.class)
    };

    @Override
    public String getParserDescription() {
        return "Provides the OU branch-rate gradient through the scaled branch-time derivative.";
    }

    @Override
    public Class getReturnType() {
        return OUBranchRateGradient.class;
    }
}
