package dr.evomodelxml.branchmodel;

import dr.evolution.tree.NodeRef;
import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import static org.apache.commons.math.stat.StatUtils.max;
import static org.apache.commons.math.stat.StatUtils.min;

/**
 * @author Filippo Monti
 */

public class RewardsAwareBranchModelParser extends AbstractXMLObjectParser {

    public final String PARSER_NAME = "rewardsAwareBranchModel";
    private final String REWARD_RATES = "rewardRates";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        SubstitutionModel underlyingSubstitutionModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);
        Parameter rewardRates = (Parameter) xo.getElementFirstChild(REWARD_RATES);

        if (rewardRates.getDimension() != underlyingSubstitutionModel.getDataType().getStateCount()) {
            throw new XMLParseException("The number of reward rates should equal to the number of states");
        }
        //     TODO   maybe this part should be moved into the RewardsAwareBranchModel constructor?

        ArbitraryBranchRates branchRateModel = (ArbitraryBranchRates) xo.getChild(BranchRateModel.class);
        TreeModel tree = (TreeModel) branchRateModel.getTree();
        checkCompatibility(rewardRates, branchRateModel, tree);

        return new RewardsAwareBranchModel(tree, underlyingSubstitutionModel, rewardRates, branchRateModel);
    }

    private void checkCompatibility(Parameter rewardRates,
                                    ArbitraryBranchRates branchRateModel,
                                    TreeModel tree) throws XMLParseException {
        double minRate = min(rewardRates.getParameterValues());
        double maxRate = max(rewardRates.getParameterValues());
        for (int i = 0; i < tree.getNodeCount(); i++) {
            final NodeRef node = tree.getNode(i);
            if (tree.isRoot(node)) continue;
            final double rate = branchRateModel.getBranchRate(tree, node);
            if (rate < minRate || rate > maxRate) {
                throw new XMLParseException("The (branch-standardized) total rewards should be within the range (min, max) of the rewards rates.");
            }
        }
    }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(BranchRateModel.class),
            new ElementRule(SubstitutionModel.class),
            new ElementRule(REWARD_RATES, Parameter.class),
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
