package dr.evomodelxml.branchmodel;

import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.RewardRates;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.branchratemodel.RewardsAwareCategoricalMixtureBranchRatesParser;
import dr.evomodelxml.branchratemodel.RewardsAwareMixtureBranchRatesParser;
import dr.inference.model.Parameter;
import dr.xml.*;


/**
 * @author Filippo Monti
 */

public class RewardsAwareBranchModelParser extends AbstractXMLObjectParser {

    public final String PARSER_NAME = "rewardsAwareBranchModel";
    private final String INDICATOR = "indicator";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        SubstitutionModel underlyingSubstitutionModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);

        ArbitraryBranchRates branchRateModel = (ArbitraryBranchRates) xo.getChild(BranchRateModel.class);
        TreeModel tree = (TreeModel) branchRateModel.getTree();

        RewardRates rewardRates = (RewardRates) xo.getChild(RewardRates.class);
        Parameter rewardRatesValues = rewardRates.getValues();
        Parameter rewardRatesMapping = rewardRates.getStateIndices();

        if (rewardRatesValues.getDimension() != underlyingSubstitutionModel.getDataType().getStateCount()) {
            throw new XMLParseException("The number of reward rates should equal to the number of states");
        }
        if (rewardRatesMapping.getDimension() != underlyingSubstitutionModel.getDataType().getStateCount()) {
            throw new XMLParseException("The reward rates mapping should have the same dimension as the number of states");
        }

        boolean conditional = xo.getAttribute( "conditional", false);

        final boolean hasLegacyState = xo.hasChildNamed(INDICATOR) ||
                xo.hasChildNamed(RewardsAwareMixtureBranchRatesParser.ATOMS_INDICES);
        final boolean hasCategoricalState = xo.hasChildNamed(
                RewardsAwareCategoricalMixtureBranchRatesParser.CATEGORY_STATE) ||
                xo.hasChildNamed(RewardsAwareCategoricalMixtureBranchRatesParser.CATEGORY_CUTS);

        if (hasLegacyState == hasCategoricalState) {
            throw new XMLParseException("Provide exactly one reward-mixture state representation: either <" +
                    INDICATOR + "> plus <" + RewardsAwareMixtureBranchRatesParser.ATOMS_INDICES + "> or <" +
                    RewardsAwareCategoricalMixtureBranchRatesParser.CATEGORY_STATE + "> plus <" +
                    RewardsAwareCategoricalMixtureBranchRatesParser.CATEGORY_CUTS + ">.");
        }

        if (hasCategoricalState) {
            if (!xo.hasChildNamed(RewardsAwareCategoricalMixtureBranchRatesParser.CATEGORY_STATE) ||
                    !xo.hasChildNamed(RewardsAwareCategoricalMixtureBranchRatesParser.CATEGORY_CUTS)) {
                throw new XMLParseException("Categorical reward-mixture state requires both <" +
                        RewardsAwareCategoricalMixtureBranchRatesParser.CATEGORY_STATE + "> and <" +
                        RewardsAwareCategoricalMixtureBranchRatesParser.CATEGORY_CUTS + ">.");
            }
            final Parameter categoryState = (Parameter) xo.getElementFirstChild(
                    RewardsAwareCategoricalMixtureBranchRatesParser.CATEGORY_STATE);
            final Parameter categoryCuts = (Parameter) xo.getElementFirstChild(
                    RewardsAwareCategoricalMixtureBranchRatesParser.CATEGORY_CUTS);

            return new RewardsAwareBranchModel(tree, underlyingSubstitutionModel,
                    rewardRates, categoryState, categoryCuts, branchRateModel, conditional);
        }

        if (!xo.hasChildNamed(INDICATOR) ||
                !xo.hasChildNamed(RewardsAwareMixtureBranchRatesParser.ATOMS_INDICES)) {
            throw new XMLParseException("Legacy reward-mixture state requires both <" +
                    INDICATOR + "> and <" + RewardsAwareMixtureBranchRatesParser.ATOMS_INDICES + ">.");
        }
        Parameter indicator = (Parameter) xo.getElementFirstChild(INDICATOR);

        Parameter atomIndices = (Parameter) xo.getElementFirstChild(RewardsAwareMixtureBranchRatesParser.ATOMS_INDICES);

        return new RewardsAwareBranchModel(tree, underlyingSubstitutionModel,
                rewardRates, indicator, branchRateModel, atomIndices, conditional);
    }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule("conditional", true),
            new ElementRule(BranchRateModel.class),
            new ElementRule(SubstitutionModel.class),
            new ElementRule(INDICATOR, Parameter.class, "Legacy 0/1 continuous-vs-atomic indicator.", true),

            new ElementRule(RewardsAwareMixtureBranchRatesParser.ATOMS_INDICES, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class)
            }, true),
            new ElementRule(RewardsAwareCategoricalMixtureBranchRatesParser.CATEGORY_STATE, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class)
            }, "Embedded categorical branch-state coordinate.", true),
            new ElementRule(RewardsAwareCategoricalMixtureBranchRatesParser.CATEGORY_CUTS, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class)
            }, "Cut points defining continuous and atomic reward categories.", true),
            new ElementRule(RewardRates.class)
    };

    @Override
    public String getParserDescription() {
        return "Parser for reward aware branch model";
    }

    @Override
    public Class getReturnType() {
        return RewardsAwareBranchModel.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

}//END: class
