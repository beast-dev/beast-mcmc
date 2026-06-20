/*
 * MixtureIndicatorRewardModelGibbsOperatorParser.java
 *
 * Parser for MixtureIndicatorRewardModelGibbsOperator
 *
 * XML:
 * <mixtureIndicatorRewardModelGibbsOperator id="..." weight="...">
 *     <indicatorZ>
 *         <parameter idref="z"/>
 *     </indicatorZ>
 *     <pi>                     (optional)
 *         <parameter idref="pi"/>
 *     </pi>
 *     <treeDataLikelihood idref="treeLikelihood"/>
 *     <rewardsAwareBranchModel idref="rewardsAwareBranchModel"/>
 *      <rewardsAwarePartialLikelihoodProvider idref="partialLikelihoodProvider"/>
 *     <updateAllBranches>true|false</updateAllBranches>   (optional, default false)
 * </mixtureIndicatorRewardModelGibbsOperator>
 */

/*
 * @author Filippo Monti
 */

package dr.inferencexml.operators;

import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptableMCMCOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.RewardsMixtureIndicatorAndAtomIndicesOperator;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

public class RewardsMixtureIndicatorAndAtomIndicesOperatorParser extends AbstractXMLObjectParser {

    public static final String OPERATOR_NAME = "rewardsMixtureIndicatorAndAtomIndicesOperator";

    private static final String INDICATOR_Z = "indicatorZ";
    private static final String ATOM_INDEX = "atomIndex";
    private static final String UPDATE_PROPORTION = "updateProportion";
    private static final String USE_CLUSTER_MOVES = "useClusterMoves";
    private static final String CLUSTER_SIZE = "clusterSize";
    private static final String CLUSTER_BORDER_BIAS = "clusterBorderBias";
    private static final String DEPENDENT_CTMC_LIKELIHOODS = "dependentCtmcLikelihoods";
    private static final String DEPENDENT_CONTINUOUS_LIKELIHOODS = "dependentContinuousLikelihoods";

    @Override
    public String getParserName() {
        return OPERATOR_NAME;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        // weight attribute is standard for operators
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final double updateProportion = xo.getDoubleAttribute(UPDATE_PROPORTION, 1);
        if (updateProportion <= 0 || updateProportion > 1) {
            throw new XMLParseException("updateProportion must be in (0,1]. Found: " + updateProportion);
        }
        final boolean adapt = xo.getBooleanAttribute(AdaptableMCMCOperator.AUTO_OPTIMIZE, false);
        final boolean useClusterMoves = xo.getBooleanAttribute(USE_CLUSTER_MOVES, false);
        final int clusterSize = xo.getIntegerAttribute(CLUSTER_SIZE, 2);
        final double clusterBorderBias = xo.getDoubleAttribute(CLUSTER_BORDER_BIAS, 0.5);

        if (clusterSize < 1) {
            throw new XMLParseException("clusterSize must be >= 1. Found: " + clusterSize);
        }
        if (clusterBorderBias < 0.0 || clusterBorderBias > 1.0) {
            throw new XMLParseException("clusterBorderBias must be in [0,1]. Found: " + clusterBorderBias);
        }

        // required: indicatorZ
        final XMLObject zxo = xo.getChild(INDICATOR_Z);
        if (zxo == null) {
            throw new XMLParseException("Missing <" + INDICATOR_Z + "> element.");
        }
        final Parameter indicatorZ = (Parameter) zxo.getChild(Parameter.class);
        final Parameter atomIndex = (Parameter) xo.getElementFirstChild(ATOM_INDEX);

        final RewardsAwareBranchModel rewardsAwareBranchModel =
                (RewardsAwareBranchModel) xo.getChild(RewardsAwareBranchModel.class);
        final TreeDataLikelihood treeDataLikelihood =
                (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
        final TreeDataLikelihood[] dependentTreeDataLikelihoods =
                parseDependentTreeDataLikelihoods(xo);
        final TreeDataLikelihood[] dependentContinuousTreeDataLikelihoods =
                parseDependentContinuousTreeDataLikelihoods(xo);

        return new RewardsMixtureIndicatorAndAtomIndicesOperator(
                indicatorZ,
                atomIndex,
                rewardsAwareBranchModel,
                treeDataLikelihood,
                dependentTreeDataLikelihoods,
                dependentContinuousTreeDataLikelihoods,
                updateProportion,
                adapt,
                useClusterMoves,
                clusterSize,
                clusterBorderBias,
                weight
        );
    }

    private TreeDataLikelihood[] parseDependentTreeDataLikelihoods(final XMLObject xo) {
        if (!xo.hasChildNamed(DEPENDENT_CTMC_LIKELIHOODS)) {
            return new TreeDataLikelihood[0];
        }

        final XMLObject dependentXo = xo.getChild(DEPENDENT_CTMC_LIKELIHOODS);
        final List<TreeDataLikelihood> likelihoods = new ArrayList<TreeDataLikelihood>();
        for (int i = 0; i < dependentXo.getChildCount(); i++) {
            final Object child = dependentXo.getChild(i);
            if (child instanceof TreeDataLikelihood) {
                likelihoods.add((TreeDataLikelihood) child);
            }
        }

        return likelihoods.toArray(new TreeDataLikelihood[likelihoods.size()]);
    }

    private TreeDataLikelihood[] parseDependentContinuousTreeDataLikelihoods(final XMLObject xo) {
        if (!xo.hasChildNamed(DEPENDENT_CONTINUOUS_LIKELIHOODS)) {
            return new TreeDataLikelihood[0];
        }

        final XMLObject dependentXo = xo.getChild(DEPENDENT_CONTINUOUS_LIKELIHOODS);
        final List<TreeDataLikelihood> likelihoods = new ArrayList<TreeDataLikelihood>();
        for (int i = 0; i < dependentXo.getChildCount(); i++) {
            final Object child = dependentXo.getChild(i);
            if (child instanceof TreeDataLikelihood) {
                likelihoods.add((TreeDataLikelihood) child);
            }
        }

        return likelihoods.toArray(new TreeDataLikelihood[likelihoods.size()]);
    }

    @Override
    public String getParserDescription() {
        return "Gibbs operator updating per-branch mixture indicator z_b selecting atomic vs continuous reward model transition matrices.";
    }

    @Override
    public Class getReturnType() {
        return RewardsMixtureIndicatorAndAtomIndicesOperator.class;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newDoubleRule(UPDATE_PROPORTION, true),
            AttributeRule.newBooleanRule(AdaptableMCMCOperator.AUTO_OPTIMIZE, true),
            AttributeRule.newBooleanRule(USE_CLUSTER_MOVES, true),
            AttributeRule.newIntegerRule(CLUSTER_SIZE, true),
            AttributeRule.newDoubleRule(CLUSTER_BORDER_BIAS, true),

            new ElementRule(RewardsAwareBranchModel.class, false),
            new ElementRule(TreeDataLikelihood.class, false),
            new ElementRule(DEPENDENT_CTMC_LIKELIHOODS,
                    new XMLSyntaxRule[] { new ElementRule(TreeDataLikelihood.class, 1, Integer.MAX_VALUE) },
                    "Optional dependent CTMC TreeDataLikelihoods sharing the reward vector.", true),
            new ElementRule(DEPENDENT_CONTINUOUS_LIKELIHOODS,
                    new XMLSyntaxRule[] { new ElementRule(TreeDataLikelihood.class, 1, Integer.MAX_VALUE) },
                    "Optional dependent continuous TreeDataLikelihoods sharing the reward vector.", true),

            new ElementRule(INDICATOR_Z,
                    new XMLSyntaxRule[] { new ElementRule(Parameter.class) },
                    "Indicator vector z_b (0/1).", false),
            new ElementRule(ATOM_INDEX,
                    new XMLSyntaxRule[] { new ElementRule(Parameter.class) },
                    "Indicator vector z_b (0/1).", false),
    };
}
