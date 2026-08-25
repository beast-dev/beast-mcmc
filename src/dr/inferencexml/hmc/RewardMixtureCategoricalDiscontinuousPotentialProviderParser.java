package dr.inferencexml.hmc;

import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.branchmodel.RewardMixtureAtomicPseudoPrior;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodelxml.branchratemodel.RewardsAwareCategoricalMixtureBranchRatesParser;
import dr.inference.hmc.DiscontinuousPotentialProvider;
import dr.inference.hmc.RewardMixtureCategoricalDiscontinuousPotentialProvider;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

import java.util.ArrayList;
import java.util.List;

/*
 * @author Filippo Monti
 */
public final class RewardMixtureCategoricalDiscontinuousPotentialProviderParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "rewardMixtureCategoricalDiscontinuousPotential";
    public static final String DEPENDENT_CTMC_LIKELIHOODS = "dependentCtmcLikelihoods";
    public static final String DEPENDENT_CONTINUOUS_LIKELIHOODS = "dependentContinuousLikelihoods";

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final Parameter categoryState = (Parameter) xo.getElementFirstChild(
                RewardsAwareCategoricalMixtureBranchRatesParser.CATEGORY_STATE);
        final Parameter categoryCuts = (Parameter) xo.getElementFirstChild(
                RewardsAwareCategoricalMixtureBranchRatesParser.CATEGORY_CUTS);
        final RewardsAwareBranchModel rewardsAwareBranchModel =
                (RewardsAwareBranchModel) xo.getChild(RewardsAwareBranchModel.class);
        final TreeDataLikelihood treeDataLikelihood =
                (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

        return new RewardMixtureCategoricalDiscontinuousPotentialProvider(
                categoryState,
                categoryCuts,
                rewardsAwareBranchModel,
                treeDataLikelihood,
                parseTreeDataLikelihoods(xo, DEPENDENT_CTMC_LIKELIHOODS),
                parseTreeDataLikelihoods(xo, DEPENDENT_CONTINUOUS_LIKELIHOODS),
                parseAtomicPseudoPrior(xo));
    }

    static RewardMixtureAtomicPseudoPrior parseAtomicPseudoPrior(final XMLObject xo) {
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
        return "Branch-local discontinuous potential for categorical reward-mixture branch states.";
    }

    @Override
    public Class getReturnType() {
        return DiscontinuousPotentialProvider.class;
    }
}
