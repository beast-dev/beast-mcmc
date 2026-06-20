package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.RewardRates;
import dr.evomodel.branchratemodel.RewardsAwareMixtureBranchRates;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRewardsAwareMixtureBranchRatesRewardRateGradient;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.hmc.NumericGradientStepSizeProvider;
import dr.xml.*;

import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;

/**
 * Parser for the continuous-trait gradient with respect to free reward-rate
 * values in a rewards-aware mixture branch-rate model.
 */
public final class ContinuousRewardsAwareMixtureBranchRatesRewardRateGradientParser
        extends AbstractXMLObjectParser {

    private static final String NAME = "continuousRewardsAwareMixtureBranchRatesRewardRateGradient";
    public static final String GRADIENT_CHECK_TOLERANCE = "gradientCheckTolerance";

    @Override
    public String getParserName() {
        return NAME;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final String traitName = xo.getAttribute(TreeTraitParserUtilities.TRAIT_NAME, DEFAULT_TRAIT_NAME);
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

        final ContinuousRewardsAwareMixtureBranchRatesRewardRateGradient gradient =
                new ContinuousRewardsAwareMixtureBranchRatesRewardRateGradient(
                        traitName,
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
            AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME, true),
            AttributeRule.newDoubleRule(GRADIENT_CHECK_TOLERANCE, true),
            AttributeRule.newStringRule(NumericGradientStepSizeProvider.NUMERIC_STEP_SIZE, true),

            new ElementRule(TreeDataLikelihood.class),
            new ElementRule(RewardsAwareMixtureBranchRates.class),
            new ElementRule(RewardRates.class)
    };

    @Override
    public String getParserDescription() {
        return "Provides the continuous-trait gradient with respect to varying reward-rate values "
                + "under RewardsAwareMixtureBranchRates.";
    }

    @Override
    public Class getReturnType() {
        return ContinuousRewardsAwareMixtureBranchRatesRewardRateGradient.class;
    }
}
