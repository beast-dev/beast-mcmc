package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.RewardsAwareMixtureBranchRates;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.OURewardsAwareMixtureBranchRatesLocationScaleGradient;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.hmc.NumericGradientStepSizeProvider;
import dr.inference.model.Parameter;
import dr.xml.*;

import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;

public final class OURewardsAwareMixtureBranchRatesLocationScaleGradientParser
        extends AbstractXMLObjectParser {

    private static final String NAME = "ouRewardsAwareMixtureBranchRatesLocationScaleGradient";
    public static final String GRADIENT_CHECK_TOLERANCE = "gradientCheckTolerance";

    @Override
    public String getParserName() {
        return NAME;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final String traitName = xo.getAttribute(TreeTraitParserUtilities.TRAIT_NAME, DEFAULT_TRAIT_NAME);
        final Double tolerance = xo.hasAttribute(GRADIENT_CHECK_TOLERANCE)
                ? xo.getDoubleAttribute(GRADIENT_CHECK_TOLERANCE)
                : null;

        final TreeDataLikelihood treeDataLikelihood =
                (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
        final RewardsAwareMixtureBranchRates branchRateModel =
                (RewardsAwareMixtureBranchRates) xo.getChild(RewardsAwareMixtureBranchRates.class);
        final Parameter locationScaleParameter = (Parameter) xo.getChild(Parameter.class);

        final OURewardsAwareMixtureBranchRatesLocationScaleGradient gradient =
                new OURewardsAwareMixtureBranchRatesLocationScaleGradient(
                        traitName,
                        treeDataLikelihood,
                        branchRateModel,
                        locationScaleParameter,
                        tolerance);

        gradient.setNumericGradientStepSize(NumericGradientStepSizeProvider.parseStepSizeRatio(xo));
        return gradient;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME, true),
            AttributeRule.newDoubleRule(GRADIENT_CHECK_TOLERANCE, true),
            AttributeRule.newStringRule(NumericGradientStepSizeProvider.NUMERIC_STEP_SIZE, true),

            new ElementRule(TreeDataLikelihood.class),
            new ElementRule(RewardsAwareMixtureBranchRates.class),
            new ElementRule(Parameter.class)
    };

    @Override
    public String getParserDescription() {
        return "Provides the OU trait gradient with respect to the location and scale "
                + "of a linear RewardsAwareMixtureBranchRates transform.";
    }

    @Override
    public Class getReturnType() {
        return OURewardsAwareMixtureBranchRatesLocationScaleGradient.class;
    }
}
