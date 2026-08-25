package dr.inferencexml.hmc;

import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.branchmodel.RewardMixtureAtomicPseudoPrior;
import dr.evomodel.branchratemodel.PerBranchRewardMixtureCategoryDecoder;
import dr.evomodel.branchratemodel.RewardMixtureCategoryDecoding;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.hmc.DiscontinuousPotentialProvider;
import dr.inference.hmc.RewardMixtureCategoricalDiscontinuousPotentialProviderDynamic;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

import java.util.ArrayList;
import java.util.List;

/**
 * Same XML surface as {@link RewardMixtureCategoricalDiscontinuousPotentialProviderParser},
 * but builds {@link RewardMixtureCategoricalDiscontinuousPotentialProviderDynamic}.
 * Unlike the static parser, this one does not read categoryState/categoryCuts
 * as separate elements: it reuses the decoder the child RewardsAwareBranchModel
 * already owns (built via rewardsAwareCategoricalMixtureBranchRatesDynamic),
 * so the two consumers share one decoder instance and one
 * refreshEmbedding() cost instead of each maintaining their own. New,
 * opt-in element name; does not modify the existing parser.
 *
 * @author Filippo Monti
 */
public final class RewardMixtureCategoricalDiscontinuousPotentialProviderDynamicParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "rewardMixtureCategoricalDiscontinuousPotentialDynamic";
    public static final String DEPENDENT_CTMC_LIKELIHOODS =
            RewardMixtureCategoricalDiscontinuousPotentialProviderParser.DEPENDENT_CTMC_LIKELIHOODS;
    public static final String DEPENDENT_CONTINUOUS_LIKELIHOODS =
            RewardMixtureCategoricalDiscontinuousPotentialProviderParser.DEPENDENT_CONTINUOUS_LIKELIHOODS;

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final RewardsAwareBranchModel rewardsAwareBranchModel =
                (RewardsAwareBranchModel) xo.getChild(RewardsAwareBranchModel.class);
        final TreeDataLikelihood treeDataLikelihood =
                (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

        final RewardMixtureCategoryDecoding decoder = rewardsAwareBranchModel.getCategoryDecoder();
        if (!(decoder instanceof PerBranchRewardMixtureCategoryDecoder)) {
            throw new XMLParseException(PARSER_NAME + " requires its <" +
                    RewardsAwareBranchModel.REWARDS_AWARE_BRANCH_MODEL +
                    "> to have been built with rewardsAwareCategoricalMixtureBranchRatesDynamic " +
                    "(a PerBranchRewardMixtureCategoryDecoder), not the shared-embedding decoder.");
        }

        return new RewardMixtureCategoricalDiscontinuousPotentialProviderDynamic(
                (PerBranchRewardMixtureCategoryDecoder) decoder,
                rewardsAwareBranchModel,
                treeDataLikelihood,
                parseTreeDataLikelihoods(xo, DEPENDENT_CTMC_LIKELIHOODS),
                parseTreeDataLikelihoods(xo, DEPENDENT_CONTINUOUS_LIKELIHOODS),
                RewardMixtureCategoricalDiscontinuousPotentialProviderParser.parseAtomicPseudoPrior(xo));
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
        return "Branch-local discontinuous potential for categorical reward-mixture branch states, with " +
                "atomic states ordered by reward value and the continuous category dynamically positioned " +
                "at each branch's live total.rewards.cts value.";
    }

    @Override
    public Class getReturnType() {
        return DiscontinuousPotentialProvider.class;
    }
}
