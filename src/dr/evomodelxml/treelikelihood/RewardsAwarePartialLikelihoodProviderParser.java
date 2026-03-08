package dr.evomodelxml.treelikelihood;

import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.RewardsAwarePartialLikelihoodProvider;
import dr.xml.*;

/*
 * Author: Filippo Monti
 */

public final class RewardsAwarePartialLikelihoodProviderParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "rewardsAwarePartialLikelihoodProvider";

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {


        TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
        RewardsAwareBranchModel rewardsAwareBranchModel = (RewardsAwareBranchModel) xo.getChild(RewardsAwareBranchModel.class);

        BeagleDataLikelihoodDelegate likelihoodDelegate = (BeagleDataLikelihoodDelegate) treeDataLikelihood.getDataLikelihoodDelegate();

        return new RewardsAwarePartialLikelihoodProvider(treeDataLikelihood, likelihoodDelegate, rewardsAwareBranchModel);
    }

    @Override
    public String getParserDescription() {
        return "Provider of pre-order partial likelihoods for rewards-aware tree likelihoods.";
    }

    @Override
    public Class getReturnType() {
        return RewardsAwarePartialLikelihoodProvider.class;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private static final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(TreeDataLikelihood.class),
            new ElementRule(RewardsAwareBranchModel.class),
    };
}
