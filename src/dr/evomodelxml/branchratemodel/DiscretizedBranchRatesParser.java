/*
 * DiscretizedBranchRatesParser.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.DiscretizedBranchRates;
import dr.evomodel.tree.TreeModel;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 */
public class DiscretizedBranchRatesParser extends AbstractXMLObjectParser {

    public static final String DISCRETIZED_BRANCH_RATES = "discretizedBranchRates";
    public static final String DISTRIBUTION = "distribution";
    public static final String RATE_CATEGORIES = "rateCategories";
    public static final String SINGLE_ROOT_RATE = "singleRootRate";
    public static final String OVERSAMPLING = "overSampling";
    public static final String NORMALIZE = "normalize";
    public static final String NORMALIZE_BRANCH_RATE_TO = "normalizeBranchRateTo";
    public static final String RANDOMIZE_RATES = "randomizeRates";
    public static final String KEEP_RATES = "keepRates";
    public static final String CACHED_RATES = "cachedRates";

    //public static final String NORMALIZED_MEAN = "normalizedMean";


    public String getParserName() {
        return DISCRETIZED_BRANCH_RATES;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final int overSampling = xo.getAttribute(OVERSAMPLING, 1);

        //final boolean normalize = xo.getBooleanAttribute(NORMALIZE, false);
        final boolean normalize = xo.getAttribute(NORMALIZE, false);
        /*if(xo.hasAttribute(NORMALIZE))
            normalize = xo.getBooleanAttribute(NORMALIZE);
        }*/
        //final double normalizeBranchRateTo = xo.getDoubleAttribute(NORMALIZE_BRANCH_RATE_TO);
        final double normalizeBranchRateTo = xo.getAttribute(NORMALIZE_BRANCH_RATE_TO, Double.NaN);

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        ParametricDistributionModel distributionModel = (ParametricDistributionModel) xo.getElementFirstChild(DISTRIBUTION);

        Parameter rateCategoryParameter = (Parameter) xo.getElementFirstChild(RATE_CATEGORIES);

        Logger.getLogger("dr.evomodel").info("\nUsing discretized relaxed clock model.");
        Logger.getLogger("dr.evomodel").info("  over sampling = " + overSampling);
        Logger.getLogger("dr.evomodel").info("  parametric model = " + distributionModel.getModelName());
        Logger.getLogger("dr.evomodel").info("   rate categories = " + rateCategoryParameter.getDimension());
        if(normalize) {
            Logger.getLogger("dr.evomodel").info("   mean rate is normalized to " + normalizeBranchRateTo);
        }

        if (xo.hasAttribute(SINGLE_ROOT_RATE)) {
            //singleRootRate = xo.getBooleanAttribute(SINGLE_ROOT_RATE);
            Logger.getLogger("dr.evomodel").warning("   WARNING: single root rate is not implemented!");
        }

        final boolean randomizeRates = xo.getAttribute(RANDOMIZE_RATES, true);
        final boolean keepRates = xo.getAttribute(KEEP_RATES, false);

        final boolean cachedRates = xo.getAttribute(CACHED_RATES, false);

        if (randomizeRates && keepRates) {
            throw new XMLParseException("Unable to both randomize and keep current rate categories");
        }

        /* if (xo.hasAttribute(NORMALIZED_MEAN)) {
            dbr.setNormalizedMean(xo.getDoubleAttribute(NORMALIZED_MEAN));
        }*/

        return new DiscretizedBranchRates(tree, rateCategoryParameter, distributionModel, overSampling, normalize,
                normalizeBranchRateTo, randomizeRates, keepRates, cachedRates);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return
                "This element returns an discretized relaxed clock model." +
                        "The branch rates are drawn from a discretized parametric distribution.";
    }

    public Class getReturnType() {
        return DiscretizedBranchRates.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(SINGLE_ROOT_RATE, true, "Whether only a single rate should be used for the two children branches of the root"),
            //AttributeRule.newDoubleRule(NORMALIZED_MEAN, true, "The mean rate to constrain branch rates to once branch lengths are taken into account"),
            AttributeRule.newIntegerRule(OVERSAMPLING, true, "The integer factor for oversampling the distribution model (1 means no oversampling)"),
            AttributeRule.newBooleanRule(NORMALIZE, true, "Whether the mean rate has to be normalized to a particular value"),
            AttributeRule.newDoubleRule(NORMALIZE_BRANCH_RATE_TO, true, "The mean rate to normalize to, if normalizing"),
            AttributeRule.newBooleanRule(RANDOMIZE_RATES, true, "Randomize initial categories"),
            AttributeRule.newBooleanRule(KEEP_RATES, true, "Keep current rate category specification"),
            AttributeRule.newBooleanRule(CACHED_RATES, true, "Cache rates between steps (default off)"),
            new ElementRule(TreeModel.class),
            new ElementRule(DISTRIBUTION, ParametricDistributionModel.class, "The distribution model for rates among branches", false),
            new ElementRule(RATE_CATEGORIES, Parameter.class, "The rate categories parameter", false),
    };
}
