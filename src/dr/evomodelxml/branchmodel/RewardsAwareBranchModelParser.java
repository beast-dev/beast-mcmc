package dr.evomodelxml.branchmodel;

import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.RewardRates;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
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

        Parameter indicator = (Parameter) xo.getElementFirstChild(INDICATOR);

        Parameter atomIndices = (Parameter) xo.getElementFirstChild(RewardsAwareMixtureBranchRatesParser.ATOMS_INDICES);
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

        return new RewardsAwareBranchModel(tree, underlyingSubstitutionModel,
                rewardRates, indicator, branchRateModel, atomIndices, conditional);
    }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule("conditional", true),
            new ElementRule(BranchRateModel.class),
            new ElementRule(SubstitutionModel.class),
            new ElementRule(INDICATOR, Parameter.class),

            new ElementRule(RewardsAwareMixtureBranchRatesParser.ATOMS_INDICES, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class)
            }),
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
