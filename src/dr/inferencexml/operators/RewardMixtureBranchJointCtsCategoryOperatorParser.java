package dr.inferencexml.operators;

import dr.evomodel.branchmodel.RewardMixtureAtomicPseudoPrior;
import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodelxml.branchratemodel.RewardsAwareCategoricalMixtureBranchRatesParser;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.RewardMixtureBranchJointCtsCategoryOperator;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * Parser for branch-local joint category/CTS reward-mixture updates.
 */
public final class RewardMixtureBranchJointCtsCategoryOperatorParser extends AbstractXMLObjectParser {

    public static final String OPERATOR_NAME = "rewardMixtureBranchJointCtsCategoryOperator";
    public static final String WINDOW_SIZE = "windowSize";
    public static final String MAX_STEPPING_OUT = "maxSteppingOut";
    public static final String MAX_SHRINK_ITERATIONS = "maxShrinkIterations";
    public static final String DEPENDENT_CTMC_LIKELIHOODS = "dependentCtmcLikelihoods";
    public static final String DEPENDENT_CONTINUOUS_LIKELIHOODS = "dependentContinuousLikelihoods";

    @Override
    public String getParserName() {
        return OPERATOR_NAME;
    }

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final double windowSize = xo.getDoubleAttribute(WINDOW_SIZE);
        final int maxSteppingOut = xo.getIntegerAttribute(MAX_STEPPING_OUT, 10);
        final int maxShrinkIterations = xo.getIntegerAttribute(MAX_SHRINK_ITERATIONS, 100);

        final Parameter categoryState = (Parameter) xo.getElementFirstChild(
                RewardsAwareCategoricalMixtureBranchRatesParser.CATEGORY_STATE);
        final Parameter categoryCuts = (Parameter) xo.getElementFirstChild(
                RewardsAwareCategoricalMixtureBranchRatesParser.CATEGORY_CUTS);
        final RewardsAwareBranchModel rewardsAwareBranchModel =
                (RewardsAwareBranchModel) xo.getChild(RewardsAwareBranchModel.class);
        final TreeDataLikelihood treeDataLikelihood =
                (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

        return new RewardMixtureBranchJointCtsCategoryOperator(
                categoryState,
                categoryCuts,
                rewardsAwareBranchModel,
                treeDataLikelihood,
                parseTreeDataLikelihoods(xo, DEPENDENT_CTMC_LIKELIHOODS),
                parseTreeDataLikelihoods(xo, DEPENDENT_CONTINUOUS_LIKELIHOODS),
                parseAtomicPseudoPrior(xo),
                windowSize,
                maxSteppingOut,
                maxShrinkIterations,
                weight);
    }

    private RewardMixtureAtomicPseudoPrior parseAtomicPseudoPrior(final XMLObject xo) {
        for (int i = 0; i < xo.getChildCount(); i++) {
            final Object child = xo.getChild(i);
            if (child instanceof RewardMixtureAtomicPseudoPrior) {
                return (RewardMixtureAtomicPseudoPrior) child;
            }
        }
        return null;
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
                AttributeRule.newDoubleRule(WINDOW_SIZE),
                AttributeRule.newIntegerRule(MAX_STEPPING_OUT, true),
                AttributeRule.newIntegerRule(MAX_SHRINK_ITERATIONS, true),
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
                new ElementRule(RewardMixtureAtomicPseudoPrior.class, 0, 1),
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
        return "Branch-local Gibbs-within-slice update for reward-mixture category and CTS coordinates.";
    }

    @Override
    public Class getReturnType() {
        return RewardMixtureBranchJointCtsCategoryOperator.class;
    }
}
