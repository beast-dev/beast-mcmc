package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.RewardsAwareMixtureBranchRates;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.RewardsAwareMixtureBranchRatesRewardRateGradient;
import dr.inference.hmc.NumericGradientStepSizeProvider;
import dr.inference.model.Parameter;
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

        final XMLObject rewardRatesObject = xo.getChild(RewardsAwareMixtureBranchRatesParser.REWARD_RATES);
        if (rewardRatesObject == null) {
            throw new XMLParseException("Must provide rewardRates.");
        }

        final Parameter rewardRatesValues = parseParameter(
                rewardRatesObject,
                RewardsAwareMixtureBranchRatesParser.REWARD_RATES_VALUES);
        final Parameter rewardRatesVaryingValues = parseParameter(
                rewardRatesObject,
                RewardsAwareMixtureBranchRatesParser.REWARD_RATES_VARYING_VALUES);

        final RewardsAwareMixtureBranchRatesRewardRateGradient gradient =
                new RewardsAwareMixtureBranchRatesRewardRateGradient(
                        treeDataLikelihood,
                        branchRateModel,
                        rewardRatesValues,
                        rewardRatesVaryingValues,
                        tolerance);

        gradient.setNumericGradientStepSize(NumericGradientStepSizeProvider.parseStepSizeRatio(xo));
        return gradient;
    }

    private static Parameter parseParameter(XMLObject xo, String name) throws XMLParseException {
        final Object child = xo.getElementFirstChild(name);
        if (!(child instanceof Parameter)) {
            throw new XMLParseException("rewardRates/" + name + " must contain a parameter.");
        }
        return (Parameter) child;
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
            new ElementRule(RewardsAwareMixtureBranchRatesParser.REWARD_RATES, new XMLSyntaxRule[]{
                    new ElementRule(RewardsAwareMixtureBranchRatesParser.REWARD_RATES_VALUES,
                            new XMLSyntaxRule[]{
                                    new ElementRule(Parameter.class)
                            }),
                    new ElementRule(RewardsAwareMixtureBranchRatesParser.REWARD_RATES_VARYING_VALUES,
                            new XMLSyntaxRule[]{
                                    new ElementRule(Parameter.class)
                            })
            })
    };

    @Override
    public String getParserDescription() {
        return "Provides the dependent-process gradient with respect to free reward-rate values "
                + "under RewardsAwareMixtureBranchRates.";
    }

    @Override
    public Class getReturnType() {
        return RewardsAwareMixtureBranchRatesRewardRateGradient.class;
    }
}
