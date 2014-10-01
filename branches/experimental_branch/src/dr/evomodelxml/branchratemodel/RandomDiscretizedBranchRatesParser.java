package dr.evomodelxml.branchratemodel;

//import dr.evomodel.branchratemodel.RandomDiscretizedBranchRates;
import dr.evomodel.tree.TreeModel;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Wai Lok Sibon Li
 */
public class RandomDiscretizedBranchRatesParser extends AbstractXMLObjectParser {

    public static final String RANDOM_DISCRETIZED_BRANCH_RATES = "randomDiscretizedBranchRates";
    public static final String DISTRIBUTION = "distribution";
    //public static final String RATE_CATEGORIES = "rateCategories";
    public static final String RATE_CATEGORY_QUANTILES = "rateCategoryQuantiles";
    public static final String SINGLE_ROOT_RATE = "singleRootRate";
    //public static final String OVERSAMPLING = "overSampling";
    public static final String NORMALIZE = "normalize";
    public static final String NORMALIZE_BRANCH_RATE_TO = "normalizeBranchRateTo";
    //public static final String NORMALIZED_MEAN = "normalizedMean";


    public String getParserName() {
        return RANDOM_DISCRETIZED_BRANCH_RATES;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

//        throw new RuntimeException("randomDiscretizedBranchRates has been renamed to continuousBranchRates (all of the " +
//                "parameters have been kept the same). Please make changes to the XML before rerunning. ");
        System.err.println("randomDiscretizedBranchRates has been renamed to continuousBranchRates (all of the " +
                        "parameters have been kept the same). Please make changes to the XML before rerunning. ");
        System.exit(1);
        return null;

        //final int overSampling = xo.getAttribute(OVERSAMPLING, 1);
//        final boolean normalize = xo.getAttribute(NORMALIZE, false);
//        final double normalizeBranchRateTo = xo.getAttribute(NORMALIZE_BRANCH_RATE_TO, Double.NaN);
//
//        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
//        ParametricDistributionModel distributionModel = (ParametricDistributionModel) xo.getElementFirstChild(DISTRIBUTION);
//
//        //Parameter rateCategoryParameter = (Parameter) xo.getElementFirstChild(RATE_CATEGORIES);
//
//        Parameter rateCategoryQuantilesParameter = (Parameter) xo.getElementFirstChild(RATE_CATEGORY_QUANTILES);
//
//        Logger.getLogger("dr.evomodel").info("Using random discretized relaxed clock model.");
//        //Logger.getLogger("dr.evomodel").info("  over sampling = " + overSampling);
//        Logger.getLogger("dr.evomodel").info("  parametric model = " + distributionModel.getModelName());
//        //Logger.getLogger("dr.evomodel").info("   rate categories = " + rateCategoryParameter.getDimension());
//        Logger.getLogger("dr.evomodel").info("   rate categories = " + rateCategoryQuantilesParameter.getDimension());
//        if(normalize) {
//            Logger.getLogger("dr.evomodel").info("   mean rate is normalized to " + normalizeBranchRateTo);
//        }
//
//        if (xo.hasAttribute(SINGLE_ROOT_RATE)) {
//            //singleRootRate = xo.getBooleanAttribute(SINGLE_ROOT_RATE);
//            Logger.getLogger("dr.evomodel").warning("   WARNING: single root rate is not implemented!");
//        }
//
//        /* if (xo.hasAttribute(NORMALIZED_MEAN)) {
//            dbr.setNormalizedMean(xo.getDoubleAttribute(NORMALIZED_MEAN));
//        }*/
//
//        return new RandomDiscretizedBranchRates(tree, /*rateCategoryParameter, */rateCategoryQuantilesParameter, distributionModel, /*overSampling,*/ normalize, normalizeBranchRateTo);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return
                "This element returns a random discretized relaxed clock model." +
                        "The branch rates are drawn from a continuous parametric distribution.";
    }

    public Class getReturnType() {
        return null; //RandomDiscretizedBranchRates.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newBooleanRule(SINGLE_ROOT_RATE, true, "Whether only a single rate should be used for the two children branches of the root"),
            //AttributeRule.newDoubleRule(NORMALIZED_MEAN, true, "The mean rate to constrain branch rates to once branch lengths are taken into account"),
            //AttributeRule.newIntegerRule(OVERSAMPLING, true, "The integer factor for oversampling the distribution model (1 means no oversampling)"),
            AttributeRule.newBooleanRule(NORMALIZE, true, "Whether the mean rate has to be normalized to a particular value"),
            AttributeRule.newDoubleRule(NORMALIZE_BRANCH_RATE_TO, true, "The mean rate to normalize to, if normalizing"),
            new ElementRule(TreeModel.class),
            new ElementRule(DISTRIBUTION, ParametricDistributionModel.class, "The distribution model for rates among branches", false),
            /*new ElementRule(RATE_CATEGORIES, Parameter.class, "The rate categories parameter", false),      */
            new ElementRule(RATE_CATEGORY_QUANTILES, Parameter.class, "The quantiles for", false),
    };
}