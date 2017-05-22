/*
 * DiscretizedLociRatesParser.java
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

package dr.oldevomodelxml.sitemodel;

import dr.xml.*;
import dr.evomodel.branchratemodel.DiscretizedBranchRates;
import dr.oldevomodel.sitemodel.DiscretizedLociRates;
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
