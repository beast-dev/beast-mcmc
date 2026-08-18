package dr.inferencexml.operators;

import dr.evomodel.branchratemodel.RewardRates;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OneZeroOneShuffleGibbsOperator;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * <oneZeroOneShuffleGibbsOperator id="..." weight="..." tol="...">
 *   <rewardRates>
 *     <fixedValues>
 *       <parameter idref="ind.rewardRates.fixedValues"/>
 *     </fixedValues>
 *     <varyingValues>
 *       <parameter idref="ind.rewardRates.varyingValues"/>
 *     </varyingValues>
 *     <stateIndices>
 *       <parameter idref="ind.rewardRates.stateIndices"/>
 *     </stateIndices>
 *   </rewardRates>
 *   <treeDataLikelihood idref="..."/>
 *   <dependentCtmcLikelihoods>              (optional)
 *     <treeDataLikelihood idref="..."/>
 *   </dependentCtmcLikelihoods>
 *   <dependentContinuousLikelihoods>        (optional)
 *     <treeDataLikelihood idref="..."/>
 *   </dependentContinuousLikelihoods>
 * </oneZeroOneShuffleGibbsOperator>
 */

/*
 * @author: Filippo Monti
 */
public final class OneZeroOneShuffleGibbsOperatorParser extends AbstractXMLObjectParser {

    public static final String OPERATOR = "oneZeroOneShuffleGibbsOperator";

    private static final String WEIGHT = "weight";
    private static final String TOL = "tol";
    private static final String DEPENDENT_CTMC_LIKELIHOODS = "dependentCtmcLikelihoods";
    private static final String DEPENDENT_CONTINUOUS_LIKELIHOODS = "dependentContinuousLikelihoods";

    @Override
    public String getParserName() {
        return OPERATOR;
    }

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {

        final double weight = xo.getDoubleAttribute(WEIGHT);
        final double tol = xo.hasAttribute(TOL) ? xo.getDoubleAttribute(TOL) : 0.0;

        final RewardRates rewardRates = (RewardRates) xo.getChild(RewardRates.class);
        final TreeDataLikelihood treeDataLikelihood =
                (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
        final TreeDataLikelihood[] dependentCtmcLikelihoods =
                parseTreeDataLikelihoods(xo, DEPENDENT_CTMC_LIKELIHOODS);
        final TreeDataLikelihood[] dependentContinuousLikelihoods =
                parseTreeDataLikelihoods(xo, DEPENDENT_CONTINUOUS_LIKELIHOODS);

        return new OneZeroOneShuffleGibbsOperator(
                rewardRates.getValues(),
                rewardRates.getStateIndices(),
                treeDataLikelihood,
                dependentCtmcLikelihoods,
                dependentContinuousLikelihoods,
                weight,
                tol);
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
    public String getParserDescription() {
        return "Full-conditional Gibbs refresh (enumerate all K! permutations) for the reward-rate state-to-slot mapping.";
    }

    @Override
    public Class getReturnType() {
        return OneZeroOneShuffleGibbsOperator.class;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private static final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(WEIGHT),
            AttributeRule.newDoubleRule(TOL, true),

            new ElementRule(RewardRates.class),
            new ElementRule(TreeDataLikelihood.class),
            new ElementRule(DEPENDENT_CTMC_LIKELIHOODS,
                    new XMLSyntaxRule[]{
                            new ElementRule(TreeDataLikelihood.class, 1, Integer.MAX_VALUE)
                    },
                    "Optional dependent CTMC TreeDataLikelihoods sharing the reward vector.", true),
            new ElementRule(DEPENDENT_CONTINUOUS_LIKELIHOODS,
                    new XMLSyntaxRule[]{
                            new ElementRule(TreeDataLikelihood.class, 1, Integer.MAX_VALUE)
                    },
                    "Optional dependent continuous TreeDataLikelihoods sharing the reward vector.", true),
    };
}
