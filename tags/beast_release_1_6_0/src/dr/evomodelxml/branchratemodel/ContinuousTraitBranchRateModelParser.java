package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.ContinuousTraitBranchRateModel;
import dr.evomodel.continuous.SampledMultivariateTraitLikelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

public class ContinuousTraitBranchRateModelParser extends AbstractXMLObjectParser {

    public static final String TRAIT_BRANCH_RATES = "continuousTraitRateModel";
    public static final String TRAIT = "trait";
    public static final String DIMENSION = "dimension";
    public static final String RATE = BranchRateModel.RATE;
    public static final String RATIO = "ratio";

    public String getParserName() {
        return TRAIT_BRANCH_RATES;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String trait = xo.getAttribute(TRAIT, "");
        int dimension = 0;
        if (xo.hasAttribute(DIMENSION)) {
            dimension = xo.getIntegerAttribute(DIMENSION) - 1;
        }

        SampledMultivariateTraitLikelihood traitLikelihood = (SampledMultivariateTraitLikelihood)
                xo.getChild(SampledMultivariateTraitLikelihood.class);
        if (traitLikelihood != null)
            trait = traitLikelihood.getTraitName();

        Logger.getLogger("dr.evomodel").info("Using trait '" + trait + "' as log rate estimates.");

        if (xo.hasChildNamed(RATE)) {
            Parameter rateParameter = (Parameter) xo.getElementFirstChild(RATE);
            Parameter ratioParameter = (Parameter) xo.getElementFirstChild(RATIO);

            return new ContinuousTraitBranchRateModel(trait, rateParameter, ratioParameter);
        } else {

            return new ContinuousTraitBranchRateModel(traitLikelihood, dimension);
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns an trait rate model." +
                        "The branch rates are an average of the rates provided by a node trait.";
    }

    public Class getReturnType() {
        return ContinuousTraitBranchRateModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
//                AttributeRule.newStringRule(TRAIT, false, "The name of the trait that provides the log rates at nodes"),
            AttributeRule.newIntegerRule(DIMENSION, true, "The dimension that supplies the rate"),
            new ElementRule(RATE, Parameter.class, "The rate parameter", true),
            new ElementRule(RATIO, Parameter.class, "The ratio parameter", true),
            new ElementRule(SampledMultivariateTraitLikelihood.class, true),
    };


}
