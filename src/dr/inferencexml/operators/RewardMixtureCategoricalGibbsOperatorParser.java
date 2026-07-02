package dr.inferencexml.operators;

import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodelxml.branchratemodel.RewardsAwareCategoricalMixtureBranchRatesParser;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.RewardMixtureCategoricalGibbsOperator;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Filippo Monti
 */
public final class RewardMixtureCategoricalGibbsOperatorParser extends AbstractXMLObjectParser {

    public static final String OPERATOR_NAME = "rewardMixtureCategoricalGibbsOperator";
    public static final String UPDATE_PROPORTION = "updateProportion";
    public static final String DEPENDENT_CTMC_LIKELIHOODS = "dependentCtmcLikelihoods";
    public static final String DEPENDENT_CONTINUOUS_LIKELIHOODS = "dependentContinuousLikelihoods";

    @Override
    public String getParserName() {
        return OPERATOR_NAME;
    }

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final double updateProportion = xo.getDoubleAttribute(UPDATE_PROPORTION, 1.0);
        if (updateProportion <= 0.0 || updateProportion > 1.0) {
            throw new XMLParseException("updateProportion must be in (0,1]. Found: " + updateProportion);
        }

        final Parameter categoryState = (Parameter) xo.getElementFirstChild(
                RewardsAwareCategoricalMixtureBranchRatesParser.CATEGORY_STATE);
        final Parameter categoryCuts = (Parameter) xo.getElementFirstChild(
                RewardsAwareCategoricalMixtureBranchRatesParser.CATEGORY_CUTS);
        final RewardsAwareBranchModel rewardsAwareBranchModel =
                (RewardsAwareBranchModel) xo.getChild(RewardsAwareBranchModel.class);
        final TreeDataLikelihood treeDataLikelihood =
                (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

        return new RewardMixtureCategoricalGibbsOperator(
                categoryState,
                categoryCuts,
                rewardsAwareBranchModel,
                treeDataLikelihood,
                parseTreeDataLikelihoods(xo, DEPENDENT_CTMC_LIKELIHOODS),
                parseTreeDataLikelihoods(xo, DEPENDENT_CONTINUOUS_LIKELIHOODS),
                updateProportion,
                weight);
    }

    private TreeDataLikelihood[] parseTreeDataLikelihoods(final XMLObject xo, final String childName) {
        if (!xo.hasChildNamed(childName)) {
            return new TreeDataLikelihood[0];
        }

        final XMLObject likelihoodsXo = xo.getChild(childName);
        final List<TreeDataLikelihood> likelihoods = new ArrayList<TreeDataLikelihood>();
        for (int i = 0; i < likelihoodsXo.getChildCount(); i++) {
            final Object child = likelihoodsXo.getChild(i);
            if (child instanceof TreeDataLikelihood) {
                likelihoods.add((TreeDataLikelihood) child);
            }
        }

        return likelihoods.toArray(new TreeDataLikelihood[likelihoods.size()]);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                AttributeRule.newDoubleRule(UPDATE_PROPORTION, true),
                new ElementRule(RewardsAwareCategoricalMixtureBranchRatesParser.CATEGORY_STATE,
                        new XMLSyntaxRule[] {
                                new ElementRule(Parameter.class)
                        }),
                new ElementRule(RewardsAwareCategoricalMixtureBranchRatesParser.CATEGORY_CUTS,
                        new XMLSyntaxRule[] {
                                new ElementRule(Parameter.class)
                        }),
                new ElementRule(RewardsAwareBranchModel.class),
                new ElementRule(TreeDataLikelihood.class),
                new ElementRule(DEPENDENT_CTMC_LIKELIHOODS,
                        new XMLSyntaxRule[] {
                                new ElementRule(TreeDataLikelihood.class, 1, Integer.MAX_VALUE)
                        },
                        "Optional dependent CTMC TreeDataLikelihoods sharing the reward vector.", true),
                new ElementRule(DEPENDENT_CONTINUOUS_LIKELIHOODS,
                        new XMLSyntaxRule[] {
                                new ElementRule(TreeDataLikelihood.class, 1, Integer.MAX_VALUE)
                        },
                        "Optional dependent continuous TreeDataLikelihoods sharing the reward vector.", true)
        };
    }

    @Override
    public String getParserDescription() {
        return "Gibbs refresh for embedded categorical reward-mixture branch states.";
    }

    @Override
    public Class getReturnType() {
        return RewardMixtureCategoricalGibbsOperator.class;
    }
}
