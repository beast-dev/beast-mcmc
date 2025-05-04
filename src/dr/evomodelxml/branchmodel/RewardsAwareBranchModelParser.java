package dr.evomodelxml.branchmodel;

import dr.evolution.tree.NodeRef;
import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import static dr.app.bss.Utils.max;
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
        double minRate = min(rewardRates.getParameterValues());
        double maxRate = max(rewardRates.getParameterValues());
        ArbitraryBranchRates branchRateModel = (ArbitraryBranchRates) xo.getChild(BranchRateModel.class);
        TreeModel tree = (TreeModel) branchRateModel.getTree();
        makeTotalRewardsCompatible(tree, branchRateModel, minRate, maxRate);

        return new RewardsAwareBranchModel(tree, underlyingSubstitutionModel, rewardRates, branchRateModel);
    }

    private void makeTotalRewardsCompatible(TreeModel tree, ArbitraryBranchRates branchRateModel,
                                            double minRate, double maxRate) {
        double minBranchRate;
        double maxBranchRate;
        double rate;

        for (int i = 0; i < branchRateModel.getTree().getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if (tree.isRoot(node)) continue;

            minBranchRate = minRate * tree.getBranchLength(node);
            maxBranchRate = maxRate * tree.getBranchLength(node);
            rate = branchRateModel.getBranchRate(tree, node);

            if (rate < minBranchRate || rate > maxBranchRate) {
                System.out.println("WARNING:Total reward *not* compatible: being set to mid point of the range");
                branchRateModel.setBranchRate(tree, node, (minBranchRate + maxBranchRate) / 2);
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
