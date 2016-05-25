/*
 * ContinuousBranchRatesParser.java
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

import dr.evomodel.branchratemodel.ContinuousBranchRates;
import dr.evomodel.tree.TreeModel;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Wai Lok Sibon Li
 */
public class ContinuousBranchRatesParser extends AbstractXMLObjectParser {

    public static final String CONTINUOUS_BRANCH_RATES = "continuousBranchRates";
    public static final String DISTRIBUTION = "distribution";
    public static final String RATE_CATEGORY_QUANTILES = "rateCategoryQuantiles";
    public static final String RATE_QUANTILES = "rateQuantiles";
    public static final String SINGLE_ROOT_RATE = "singleRootRate";
    public static final String NORMALIZE = "normalize";
    public static final String NORMALIZE_BRANCH_RATE_TO = "normalizeBranchRateTo";


    public String getParserName() {
        return CONTINUOUS_BRANCH_RATES;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final boolean normalize = xo.getAttribute(NORMALIZE, false);
        final double normalizeBranchRateTo = xo.getAttribute(NORMALIZE_BRANCH_RATE_TO, Double.NaN);

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        ParametricDistributionModel distributionModel = (ParametricDistributionModel) xo.getElementFirstChild(DISTRIBUTION);

        Parameter rateQuantilesParameter;
        if (xo.hasChildNamed(RATE_QUANTILES)) {
            rateQuantilesParameter = (Parameter) xo.getElementFirstChild(RATE_QUANTILES);
        } else {
            rateQuantilesParameter = (Parameter) xo.getElementFirstChild(RATE_CATEGORY_QUANTILES);
        }

        Logger.getLogger("dr.evomodel").info("Using continuous relaxed clock model.");
        Logger.getLogger("dr.evomodel").info("  parametric model = " + distributionModel.getModelName());
        Logger.getLogger("dr.evomodel").info("   rate categories = " + rateQuantilesParameter.getDimension());
        if(normalize) {
            Logger.getLogger("dr.evomodel").info("   mean rate is normalized to " + normalizeBranchRateTo);
        }

        if (xo.hasAttribute(SINGLE_ROOT_RATE)) {
            Logger.getLogger("dr.evomodel").warning("   WARNING: single root rate is not implemented!");
        }

        return new ContinuousBranchRates(tree, rateQuantilesParameter, distributionModel, normalize, normalizeBranchRateTo);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return
                "This element returns a continuous quantile uncorrelated relaxed clock model.";
    }

    public Class getReturnType() {
        return ContinuousBranchRates.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newBooleanRule(SINGLE_ROOT_RATE, true, "Whether only a single rate should be used for the two children branches of the root"),
            AttributeRule.newBooleanRule(NORMALIZE, true, "Whether the mean rate has to be normalized to a particular value"),
            AttributeRule.newDoubleRule(NORMALIZE_BRANCH_RATE_TO, true, "The mean rate to normalize to, if normalizing"),
            new ElementRule(TreeModel.class),
            new ElementRule(DISTRIBUTION, ParametricDistributionModel.class, "The distribution model for rates among branches", false),
            new XORRule(
                    new ElementRule(RATE_QUANTILES, Parameter.class, "The quantiles for each branch", false),
                    new ElementRule(RATE_CATEGORY_QUANTILES, Parameter.class, "The quantiles for each branch", false)
            )
    };
}