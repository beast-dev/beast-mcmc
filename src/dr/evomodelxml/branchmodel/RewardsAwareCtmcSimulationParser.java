package dr.evomodelxml.branchmodel;

import dr.evomodel.branchmodel.RewardsAwareCtmcSimulation;
import dr.evomodel.branchratemodel.RewardRates;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.branchratemodel.RewardsAwareMixtureBranchRatesParser;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;
import dr.xml.XORRule;

/*
 * @author Filippo Monti
 */
public final class RewardsAwareCtmcSimulationParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "rewardsAwareCtmcSimulation";
    public static final String STATE_COUNT = "stateCount";
    public static final String ROOT_STATE = "rootState";
    public static final String SEED = "seed";
    public static final String Q_MATRIX = "qMatrix";
    public static final String REWARD_RATES = "rewardRates";
    public static final String ROOT_FREQUENCIES = "rootFrequencies";
    public static final String BRANCH_REWARD_TOTALS = "branchRewardTotals";

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
        if (xo.hasAttribute(SEED)) {
            MathUtils.setSeed(xo.getIntegerAttribute(SEED));
        }

        final TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        final Parameter qMatrixParameter = parseParameter(xo, Q_MATRIX, true);
        final Parameter rewardRatesParameter = parseRewardRatesParameter(xo);
        final double[] qMatrix = qMatrixParameter.getParameterValues();
        final double[] rewardRates = rewardRatesParameter.getParameterValues();

        final int stateCount = xo.hasAttribute(STATE_COUNT)
                ? xo.getIntegerAttribute(STATE_COUNT)
                : rewardRates.length;

        final boolean hasRootState = xo.hasAttribute(ROOT_STATE);
        final boolean hasRootFrequencies = xo.hasChildNamed(ROOT_FREQUENCIES);
        if (hasRootState && hasRootFrequencies) {
            throw new XMLParseException(PARSER_NAME + " cannot specify both rootState and rootFrequencies.");
        }

        final RewardsAwareCtmcSimulation simulation;
        if (hasRootState) {
            simulation = RewardsAwareCtmcSimulation.withFixedRootState(
                    xo.hasId() ? xo.getId() : PARSER_NAME,
                    stateCount,
                    rewardRates,
                    tree,
                    qMatrix,
                    xo.getIntegerAttribute(ROOT_STATE)
            );
        } else if (hasRootFrequencies) {
            final Parameter rootFrequencies = parseParameter(xo, ROOT_FREQUENCIES, true);
            simulation = RewardsAwareCtmcSimulation.withRootFrequencies(
                    xo.hasId() ? xo.getId() : PARSER_NAME,
                    stateCount,
                    rewardRates,
                    tree,
                    qMatrix,
                    rootFrequencies.getParameterValues()
            );
        } else {
            simulation = new RewardsAwareCtmcSimulation(
                    xo.hasId() ? xo.getId() : PARSER_NAME,
                    stateCount,
                    rewardRates,
                    tree,
                    qMatrix,
                    null,
                    -1
            );
        }

        final Parameter branchRewardTotals = parseParameter(xo, BRANCH_REWARD_TOTALS, false);
        final Parameter branchRewardProportions =
                parseParameter(xo, RewardsAwareMixtureBranchRatesParser.CTS, false);
        final Parameter indicators =
                parseParameter(xo, RewardsAwareMixtureBranchRatesParser.INDICATOR, false);
        final Parameter atomIndices =
                parseParameter(xo, RewardsAwareMixtureBranchRatesParser.ATOMS_INDICES, false);

        simulation.copyBranchRewardTotalsInto(branchRewardTotals);
        simulation.copyBranchRewardProportionsInto(branchRewardProportions);
        simulation.copyIndicatorsInto(indicators);
        simulation.copyAtomIndicesInto(atomIndices);

        return simulation;
    }

    private static Parameter parseRewardRatesParameter(final XMLObject xo) throws XMLParseException {
        final RewardRates rewardRates = (RewardRates) xo.getChild(RewardRates.class);
        if (rewardRates != null) {
            return rewardRates.getValues();
        }
        return parseParameter(xo, REWARD_RATES, true);
    }

    private static Parameter parseParameter(final XMLObject xo,
                                            final String elementName,
                                            final boolean required) throws XMLParseException {
        if (!xo.hasChildNamed(elementName)) {
            if (required) {
                throw new XMLParseException(PARSER_NAME + " requires a " + elementName + " element.");
            }
            return null;
        }
        final Object child = xo.getElementFirstChild(elementName);
        if (!(child instanceof Parameter)) {
            throw new XMLParseException(PARSER_NAME + "/" + elementName + " must contain a parameter.");
        }
        return (Parameter) child;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private static final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(STATE_COUNT, true),
            AttributeRule.newIntegerRule(ROOT_STATE, true),
            AttributeRule.newIntegerRule(SEED, true),
            new ElementRule(TreeModel.class),
            new ElementRule(Q_MATRIX, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new XORRule(
                    new ElementRule(REWARD_RATES, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }),
                    new ElementRule(RewardRates.class)
            ),
            new ElementRule(ROOT_FREQUENCIES, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
            new ElementRule(BRANCH_REWARD_TOTALS, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
            new ElementRule(RewardsAwareMixtureBranchRatesParser.CTS, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
            new ElementRule(RewardsAwareMixtureBranchRatesParser.INDICATOR, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
            new ElementRule(RewardsAwareMixtureBranchRatesParser.ATOMS_INDICES, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true)
    };

    @Override
    public String getParserDescription() {
        return "Simulates one realized CTMC history on a tree and exports reward-aware branch rewards, indicators, atom indices, generated tips, and jump summaries.";
    }

    @Override
    public Class getReturnType() {
        return RewardsAwareCtmcSimulation.class;
    }
}
