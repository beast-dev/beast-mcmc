package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.RewardRates;
import dr.evomodel.branchratemodel.RewardsAwareCategoricalMixtureBranchRatesDynamic;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Same XML surface as {@link RewardsAwareCategoricalMixtureBranchRatesParser},
 * but builds {@link RewardsAwareCategoricalMixtureBranchRatesDynamic}: atomic
 * states are ordered along the embedded axis by reward value and the
 * continuous category tracks each branch's live total.rewards.cts value,
 * instead of one embedding shared by every branch (see
 * REWARD_CATEGORY_DYNAMIC_ORDERING_PLAN.md). New, opt-in element name; does
 * not modify the existing parser.
 *
 * @author Filippo Monti
 */
public final class RewardsAwareCategoricalMixtureBranchRatesDynamicParser extends AbstractXMLObjectParser {

    public static final String NAME = "rewardsAwareCategoricalMixtureBranchRatesDynamic";
    public static final String CTS = RewardsAwareCategoricalMixtureBranchRatesParser.CTS;
    public static final String CATEGORY_STATE = RewardsAwareCategoricalMixtureBranchRatesParser.CATEGORY_STATE;
    public static final String CATEGORY_CUTS = RewardsAwareCategoricalMixtureBranchRatesParser.CATEGORY_CUTS;
    public static final String INCLUDE_ROOT = RewardsAwareCategoricalMixtureBranchRatesParser.INCLUDE_ROOT;

    @Override
    public String getParserName() {
        return NAME;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        final Parameter cts = (Parameter) xo.getElementFirstChild(CTS);
        final Parameter categoryState = (Parameter) xo.getElementFirstChild(CATEGORY_STATE);
        final Parameter categoryCuts = (Parameter) xo.getElementFirstChild(CATEGORY_CUTS);
        final RewardRates rewardRates = (RewardRates) xo.getChild(RewardRates.class);

        final TreeParameterModel.Type includeRoot = xo.getAttribute(INCLUDE_ROOT, false)
                ? TreeParameterModel.Type.WITH_ROOT
                : TreeParameterModel.Type.WITHOUT_ROOT;

        final int branchCount = tree.getNodeCount() - 1;
        final int expectedDimension = includeRoot == TreeParameterModel.Type.WITH_ROOT
                ? branchCount + 1
                : branchCount;

        if (cts.getDimension() != expectedDimension) {
            throw new XMLParseException("cts dim must be " +
                    (includeRoot == TreeParameterModel.Type.WITH_ROOT ? "nodeCount" : "nodeCount-1"));
        }
        if (categoryState.getDimension() != expectedDimension) {
            throw new XMLParseException("categoryState dim must be " +
                    (includeRoot == TreeParameterModel.Type.WITH_ROOT ? "nodeCount" : "nodeCount-1"));
        }

        final int expectedCuts = rewardRates.getStateIndices().getDimension() + 2;
        if (categoryCuts.getDimension() != expectedCuts) {
            throw new XMLParseException("categoryCuts dim must be stateCount + 2 (" + expectedCuts + ")");
        }

        final ArbitraryBranchRates.BranchRateTransform transform =
                ArbitraryBranchRatesParser.parseTransform(xo);

        return new RewardsAwareCategoricalMixtureBranchRatesDynamic(
                tree,
                cts,
                categoryState,
                categoryCuts,
                rewardRates,
                transform,
                false,
                includeRoot);
    }

    @Override
    public String getParserDescription() {
        return "Reward-mixture branch-rate model driven by one embedded categorical coordinate per branch, " +
                "with atomic states ordered by reward value and the continuous category dynamically " +
                "positioned at each branch's live cts value.";
    }

    @Override
    public Class getReturnType() {
        return RewardsAwareCategoricalMixtureBranchRatesDynamic.class;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                new ElementRule(TreeModel.class),
                new ElementRule(CTS, new XMLSyntaxRule[] {
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(CATEGORY_STATE, new XMLSyntaxRule[] {
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(CATEGORY_CUTS, new XMLSyntaxRule[] {
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(RewardRates.class),
                AttributeRule.newBooleanRule(INCLUDE_ROOT, true),
                AttributeRule.newBooleanRule(ArbitraryBranchRatesParser.LINEAR_LOCATION_SCALE, true),
                new ElementRule(ArbitraryBranchRatesParser.SCALE,
                        Parameter.class, "optional scale parameter", true),
                new ElementRule(ArbitraryBranchRatesParser.LOCATION,
                        Parameter.class, "optional location parameter", true),
                new ElementRule(ArbitraryBranchRatesParser.RANDOM_INDICATOR,
                        new XMLSyntaxRule[] {
                                new ElementRule(Parameter.class),
                        }, true),
                AttributeRule.newBooleanRule(ArbitraryBranchRatesParser.SHRINKAGE, true),
        };
    }
}
