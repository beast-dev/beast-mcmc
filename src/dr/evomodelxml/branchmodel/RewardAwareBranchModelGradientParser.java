package dr.evomodelxml.branchmodel;

import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.branchmodel.RewardsAwareBranchModelGradient;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.RewardsAwarePartialLikelihoodProvider;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.Parameter;
import dr.xml.*;

import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;

/**
 * XML sketch:
 *
 * <rewardAwareBranchModelGradient traitName="rewardGrad">
 *     <treeDataLikelihood idref="treeLikelihood"/>
 *     <rewardsAwareBranchModel idref="rewardsBranchModel"/>
 *     <arbitraryBranchRates idref="ind.totalRewardsBranchRateModel"/>
 *     <rewardsAwarePartialLikelihoodProvider id="partialLikelihoodProvider"/>
 * </rewardAwareBranchModelGradient>
 *
 * @author Filippo Monti
 */
public class RewardAwareBranchModelGradientParser extends AbstractXMLObjectParser {

    private static final String NAME = "rewardAwareBranchModelGradient";

    private static final String TRAIT_NAME = TreeTraitParserUtilities.TRAIT_NAME;
    public static final String GRADIENT_CHECK_TOLERANCE = "gradientCheckTolerance";
    public static final String INDICATOR = "indicator";
    public static final String IS_INDICATORBASE_ONE = "isIndicatorBaseOne";

    @Override
    public String getParserName() {
        return NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final String traitName = xo.getAttribute(TRAIT_NAME, DEFAULT_TRAIT_NAME);
        final boolean indicatorBaseOneStr = xo.getBooleanAttribute(IS_INDICATORBASE_ONE, true);

        final Double tolerance = xo.hasAttribute(GRADIENT_CHECK_TOLERANCE)
                ? xo.getDoubleAttribute(GRADIENT_CHECK_TOLERANCE)
                : null;

        final TreeDataLikelihood treeDataLikelihood =
                (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

        final RewardsAwareBranchModel rewardsAwareBranchModel =
                (RewardsAwareBranchModel) xo.getChild(RewardsAwareBranchModel.class);

        final ArbitraryBranchRates totalRewardsBranchRates =
                (ArbitraryBranchRates) xo.getChild(ArbitraryBranchRates.class);


        final Parameter indicator = (Parameter) xo.getElementFirstChild(INDICATOR);

            return new RewardsAwareBranchModelGradient(
                    treeDataLikelihood,
                    rewardsAwareBranchModel,
                    totalRewardsBranchRates,
                    indicator,
                    tolerance,
                    false, // useHessian  TODO add this
                    indicatorBaseOneStr
            );

    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(TRAIT_NAME, true),
            AttributeRule.newDoubleRule(GRADIENT_CHECK_TOLERANCE, true),
            AttributeRule.newBooleanRule(IS_INDICATORBASE_ONE, true),

            new ElementRule(TreeDataLikelihood.class),
            new ElementRule(RewardsAwareBranchModel.class),
            new ElementRule(ArbitraryBranchRates.class),
            new ElementRule(INDICATOR, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class),
            }),
    };


    @Override
    public String getParserDescription() {
        return "Provides gradient of the TreeDataLikelihood with respect to branch-specific total rewards "
                + "(stored in ArbitraryBranchRates) under RewardsAwareBranchModel.";
    }

    @Override
    public Class getReturnType() {
        return RewardsAwareBranchModelGradient.class;
    }
}
