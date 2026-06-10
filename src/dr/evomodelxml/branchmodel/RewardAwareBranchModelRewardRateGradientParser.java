package dr.evomodelxml.branchmodel;

import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.branchmodel.RewardsAwareBranchModelRewardRateGradient;
import dr.evomodel.branchratemodel.RewardRates;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.hmc.NumericGradientStepSizeProvider;
import dr.xml.*;

public class RewardAwareBranchModelRewardRateGradientParser extends AbstractXMLObjectParser {

    private static final String NAME = "rewardAwareBranchModelRewardRateGradient";
    public static final String GRADIENT_CHECK_TOLERANCE = "gradientCheckTolerance";

    @Override
    public String getParserName() {
        return NAME;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final Double tolerance = xo.hasAttribute(GRADIENT_CHECK_TOLERANCE)
                ? xo.getDoubleAttribute(GRADIENT_CHECK_TOLERANCE)
                : null;

        final TreeDataLikelihood treeDataLikelihood =
                (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
        final RewardsAwareBranchModel rewardsAwareBranchModel =
                (RewardsAwareBranchModel) xo.getChild(RewardsAwareBranchModel.class);

        final RewardRates rewardRates = (RewardRates) xo.getChild(RewardRates.class);
        if (rewardRates.getValues() != rewardsAwareBranchModel.getRewardRatesValues() ||
                rewardRates.getVaryingValues() != rewardsAwareBranchModel.getRewardRatesInternal()) {
            throw new XMLParseException("rewardRates must reference the reward-rate parameters used by rewardsAwareBranchModel.");
        }

        final RewardsAwareBranchModelRewardRateGradient gradient =
                new RewardsAwareBranchModelRewardRateGradient(
                        treeDataLikelihood,
                        rewardsAwareBranchModel,
                        rewardRates.getValues(),
                        rewardRates.getVaryingValues(),
                        tolerance);

        gradient.setNumericGradientStepSize(NumericGradientStepSizeProvider.parseStepSizeRatio(xo));
        return gradient;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(GRADIENT_CHECK_TOLERANCE, true),
            AttributeRule.newStringRule(NumericGradientStepSizeProvider.NUMERIC_STEP_SIZE, true),

            new ElementRule(TreeDataLikelihood.class),
            new ElementRule(RewardsAwareBranchModel.class),
            new ElementRule(RewardRates.class)
    };

    @Override
    public String getParserDescription() {
        return "Provides gradient of the TreeDataLikelihood with respect to varying reward-rate values "
                + "under RewardsAwareBranchModel.";
    }

    @Override
    public Class getReturnType() {
        return RewardsAwareBranchModelRewardRateGradient.class;
    }
}
