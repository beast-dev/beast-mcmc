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
import dr.evomodel.treedatalikelihood.preorder.RewardsAwarePartialLikelihoodProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptableMCMCOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.RewardsMixtureIndicatorAndAtomIndicesOperator;
import dr.xml.*;

public class RewardsMixtureIndicatorAndAtomIndicesOperatorParser extends AbstractXMLObjectParser {

    public static final String OPERATOR_NAME = "rewardsMixtureIndicatorAndAtomIndicesOperator";

    private static final String INDICATOR_Z = "indicatorZ";
    private static final String ATOM_INDEX = "atomIndex";
    private static final String PI = "pi";
    private static final String UPDATE_PROPORTION = "updateProportion";

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

        // required: indicatorZ
        final XMLObject zxo = xo.getChild(INDICATOR_Z);
        if (zxo == null) {
            throw new XMLParseException("Missing <" + INDICATOR_Z + "> element.");
        }
        final Parameter indicatorZ = (Parameter) zxo.getChild(Parameter.class);
        final Parameter atomIndex = (Parameter) xo.getElementFirstChild(ATOM_INDEX);

        // optional: pi
        Parameter pi = null;
        final XMLObject pix = xo.getChild(PI);
        if (pix != null) {
            pi = (Parameter) pix.getChild(Parameter.class);
        }

        final RewardsAwareBranchModel rewardsAwareBranchModel =
                (RewardsAwareBranchModel) xo.getChild(RewardsAwareBranchModel.class);
        final RewardsAwarePartialLikelihoodProvider partialLikelihoodProvider =
                (RewardsAwarePartialLikelihoodProvider) xo.getChild(RewardsAwarePartialLikelihoodProvider.class);

        return new RewardsMixtureIndicatorAndAtomIndicesOperator(
                indicatorZ,
                atomIndex,
                pi,
                rewardsAwareBranchModel,
                partialLikelihoodProvider,
                updateProportion,
                adapt,
                weight
        );
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

            new ElementRule(RewardsAwarePartialLikelihoodProvider.class, false),
            new ElementRule(RewardsAwareBranchModel.class, false),

            new ElementRule(INDICATOR_Z,
                    new XMLSyntaxRule[] { new ElementRule(Parameter.class) },
                    "Indicator vector z_b (0/1).", false),
            new ElementRule(ATOM_INDEX,
                    new XMLSyntaxRule[] { new ElementRule(Parameter.class) },
                    "Indicator vector z_b (0/1).", false),
            new ElementRule(PI,
                    new XMLSyntaxRule[] { new ElementRule(Parameter.class) },
                    "Optional weights pi over states (dimension K). If omitted, uniform weights are used.",
                    true),
    };
}