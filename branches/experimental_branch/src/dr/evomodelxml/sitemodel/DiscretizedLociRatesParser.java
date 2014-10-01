package dr.evomodelxml.sitemodel;

import dr.xml.*;
import dr.evomodel.branchratemodel.DiscretizedBranchRates;
import dr.evomodel.sitemodel.DiscretizedLociRates;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Parameter;
import dr.inference.model.CompoundParameter;

import java.util.logging.Logger;

/**
 * @author Chieh-Hsi Wu
 * Parser of DiscretizedLociRates
 */
public class DiscretizedLociRatesParser extends AbstractXMLObjectParser {
    public final String DISCRETIZED_LOCI_RATES = "DiscretizedLociRates";
    public final String NORMALIZE  = "normalize";
    public final String NORMALIZE_MEAN_LOCI_RATE_TO = "normalizeMeanLociRateTo";
    public final String RATE_CATEGORIES = "rateCategories";
    public final String CATEGORY_COUNT = "categoryCount";
    public final String LOCI_RATES = "lociRates";
    public final String DISTRIBUTION = "distribution";



    public String getParserName(){
        return DISCRETIZED_LOCI_RATES;
    }
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {


        final boolean normalize = xo.getAttribute(NORMALIZE, false);
        final double normalizeBranchRateTo = xo.getAttribute(NORMALIZE_MEAN_LOCI_RATE_TO, Double.NaN);
        final int categoryCount = xo.getIntegerAttribute(CATEGORY_COUNT);

        CompoundParameter lociRates = (CompoundParameter) xo.getElementFirstChild(LOCI_RATES);
        Parameter rateCategoryParameter = (Parameter) xo.getElementFirstChild(RATE_CATEGORIES);
        ParametricDistributionModel distributionModel = (ParametricDistributionModel) xo.getElementFirstChild(DISTRIBUTION);



        Logger.getLogger("dr.evomodel").info("Using discretized loci rates model.");
        Logger.getLogger("dr.evomodel").info("Number of categories: "+categoryCount);
        Logger.getLogger("dr.evomodel").info("  parametric model = " + distributionModel.getModelName());
        if(normalize) {
            Logger.getLogger("dr.evomodel").info("   mean rate is normalized to " + normalizeBranchRateTo);
        }


        return new DiscretizedLociRates(
                lociRates,
                rateCategoryParameter,
                distributionModel,
                normalize,
                normalizeBranchRateTo,
                categoryCount);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return
                "This element returns an discretized loci rate model." +
                        "The loci rates are drawn from a discretized parametric distribution.";
    }

    public Class getReturnType() {
        return DiscretizedBranchRates.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(NORMALIZE, true, "Whether the mean rate has to be normalized to a particular value"),
            AttributeRule.newDoubleRule(NORMALIZE_MEAN_LOCI_RATE_TO, true, "The mean rate to normalize to, if normalizing"),
            AttributeRule.newIntegerRule(CATEGORY_COUNT, true, "The number of categories that the distribution will be divided into"),
            new ElementRule(LOCI_RATES, CompoundParameter.class, "The compound parameter that contains all the loci rate parameters", false),
            new ElementRule(DISTRIBUTION, ParametricDistributionModel.class, "The distribution model for rates among branches", false),
            new ElementRule(RATE_CATEGORIES, Parameter.class, "The rate categories parameter", false),
    };
}
