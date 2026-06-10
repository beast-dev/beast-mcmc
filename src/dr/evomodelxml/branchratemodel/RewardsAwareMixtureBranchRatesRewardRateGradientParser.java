package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.RewardRates;
import dr.evomodel.branchratemodel.RewardsAwareMixtureBranchRates;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.RewardsAwareMixtureBranchRatesRewardRateGradient;
import dr.inference.hmc.NumericGradientStepSizeProvider;
import dr.xml.*;

/**
 * Parser for the dependent-process gradient with respect to free reward-rate
 * values in a rewards-aware mixture branch-rate model.
 *
 * @author Filippo Monti
 */
public final class RewardsAwareMixtureBranchRatesRewardRateGradientParser extends AbstractXMLObjectParser {

    private static final String NAME = "rewardsAwareMixtureBranchRatesRewardRateGradient";
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
        final RewardsAwareMixtureBranchRates branchRateModel =
                (RewardsAwareMixtureBranchRates) xo.getChild(RewardsAwareMixtureBranchRates.class);

        final RewardRates rewardRates = (RewardRates) xo.getChild(RewardRates.class);
        if (rewardRates.getValues() != branchRateModel.getRewardRatesValues() ||
                rewardRates.getVaryingValues() != branchRateModel.getRewardRatesInternal() ||
                rewardRates.getStateIndices() != branchRateModel.getRewardRatesMapping()) {
            throw new XMLParseException("rewardRates must reference the reward-rate parameters used by rewardsAwareMixtureBranchRates.");
        }

        final RewardsAwareMixtureBranchRatesRewardRateGradient gradient =
                new RewardsAwareMixtureBranchRatesRewardRateGradient(
                        treeDataLikelihood,
                        branchRateModel,
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
            new ElementRule(RewardsAwareMixtureBranchRates.class),
            new ElementRule(RewardRates.class)
    };

    @Override
    public String getParserDescription() {
        return "Provides the dependent-process gradient with respect to varying reward-rate values "
                + "under RewardsAwareMixtureBranchRates.";
    }

    @Override
    public Class getReturnType() {
        return RewardsAwareMixtureBranchRatesRewardRateGradient.class;
    }
}
