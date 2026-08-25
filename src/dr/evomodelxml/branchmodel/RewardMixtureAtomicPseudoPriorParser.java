package dr.evomodelxml.branchmodel;

import dr.evomodel.branchmodel.RewardMixtureAtomicPseudoPrior;
import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.inference.hmc.NumericGradientStepSizeProvider;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Parser for the reward-mixture atomic pseudo-prior.
 *
 * @author Filippo Monti
 */
public final class RewardMixtureAtomicPseudoPriorParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "rewardMixtureAtomicPseudoPrior";
    public static final String STANDARD_DEVIATION = "standardDeviation";

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
        final RewardsAwareBranchModel rewardsAwareBranchModel =
                (RewardsAwareBranchModel) xo.getChild(RewardsAwareBranchModel.class);
        final ArbitraryBranchRates totalRewardsBranchRates =
                (ArbitraryBranchRates) xo.getChild(ArbitraryBranchRates.class);
        final double standardDeviation = xo.getDoubleAttribute(STANDARD_DEVIATION);

        final RewardMixtureAtomicPseudoPrior pseudoPrior =
                new RewardMixtureAtomicPseudoPrior(
                        rewardsAwareBranchModel,
                        totalRewardsBranchRates,
                        standardDeviation);
        pseudoPrior.setNumericGradientStepSize(NumericGradientStepSizeProvider.parseStepSizeRatio(xo));
        return pseudoPrior;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(STANDARD_DEVIATION),
            AttributeRule.newStringRule(NumericGradientStepSizeProvider.NUMERIC_STEP_SIZE, true),
            new ElementRule(RewardsAwareBranchModel.class),
            new ElementRule(ArbitraryBranchRates.class)
    };

    @Override
    public String getParserDescription() {
        return "Conditional pseudo-prior for total.rewards.cts in atomic reward-mixture categories.";
    }

    @Override
    public Class getReturnType() {
        return RewardMixtureAtomicPseudoPrior.class;
    }
}
