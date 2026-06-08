package dr.evomodelxml.branchmodel;

import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.branchmodel.RewardsAwareBranchModelRewardRateGradient;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodelxml.branchratemodel.RewardsAwareMixtureBranchRatesParser;
import dr.inference.hmc.NumericGradientStepSizeProvider;
import dr.inference.model.Parameter;
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

        final Parameter rewardRatesValues = parseRewardRatesValues(xo);

        final RewardsAwareBranchModelRewardRateGradient gradient =
                new RewardsAwareBranchModelRewardRateGradient(
                        treeDataLikelihood,
                        rewardsAwareBranchModel,
                        rewardRatesValues,
                        tolerance);

        gradient.setNumericGradientStepSize(NumericGradientStepSizeProvider.parseStepSizeRatio(xo));
        return gradient;
    }

    private Parameter parseRewardRatesValues(XMLObject xo) throws XMLParseException {
        final XMLObject rewardRatesObject = xo.getChild(RewardsAwareMixtureBranchRatesParser.REWARD_RATES);
        if (rewardRatesObject != null) {
            final Object child = rewardRatesObject.getElementFirstChild(
                    RewardsAwareMixtureBranchRatesParser.REWARD_RATES_VALUES);
            if (!(child instanceof Parameter)) {
                throw new XMLParseException("rewardRates/values must contain a parameter.");
            }
            return (Parameter) child;
        }

        throw new XMLParseException("Must provide rewardRates/values.");
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
            new ElementRule(RewardsAwareMixtureBranchRatesParser.REWARD_RATES, new XMLSyntaxRule[]{
                    new ElementRule(RewardsAwareMixtureBranchRatesParser.REWARD_RATES_VALUES,
                            new XMLSyntaxRule[]{
                                    new ElementRule(Parameter.class)
                            })
            })
    };

    @Override
    public String getParserDescription() {
        return "Provides gradient of the TreeDataLikelihood with respect to reward-rate values "
                + "under RewardsAwareBranchModel.";
    }

    @Override
    public Class getReturnType() {
        return RewardsAwareBranchModelRewardRateGradient.class;
    }
}
