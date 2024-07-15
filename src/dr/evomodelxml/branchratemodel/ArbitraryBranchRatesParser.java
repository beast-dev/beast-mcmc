/*
 * ArbitraryBranchRatesParser.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchSpecificFixedEffects;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.distribution.shrinkage.BayesianBridgeLikelihood;
import dr.inference.distribution.shrinkage.BayesianBridgeStatisticsProvider;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.logging.Logger;

import static dr.evomodel.branchratemodel.ArbitraryBranchRates.make;

/**
 */
public class ArbitraryBranchRatesParser extends AbstractXMLObjectParser {

    public static final String ARBITRARY_BRANCH_RATES = "arbitraryBranchRates";
    public static final String RATES = "rates";
    public static final String RECIPROCAL = "reciprocal";
    public static final String EXP = "exp";
    public static final String MULTIPLIER = "multiplier";
    public static final String CENTER_AT_ONE = "centerAtOne";
    public static final String RANDOMIZE_RATES = "randomizeRates";
    public  static final String RANDOM_SCALE = "randomScale";

    public static final String INCLUDE_ROOT = "includeRoot";
    public static final String RANDOM_INDICATOR = "randomIndicator"; // keep some rates fixed but randomize others

    public static final String SHRINKAGE = "shrinkage";

    public static final String LOCATION = "location";
    public static final String SCALE = "scale";

    public String getParserName() {
        return ARBITRARY_BRANCH_RATES;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

        XMLObject cxo = xo.getChild(RATES);

        Parameter rateCategoryParameter = (Parameter) cxo.getChild(Parameter.class);

        boolean centerAtOne = xo.getAttribute(CENTER_AT_ONE, true);

        boolean randomizeRates = xo.getAttribute(RANDOMIZE_RATES, false);

//        if (centerAtOne && randomizeRates) {
//            throw new XMLParseException("Cannot centerAtOne and randomize the starting rates");
//        }

        final int numBranches = tree.getNodeCount() - 1;
        if (rateCategoryParameter.getDimension() > 1 && (rateCategoryParameter.getDimension() != numBranches)) {
            throw new XMLParseException("Incorrect number of rate parameters");
        }

        if (rateCategoryParameter.getDimension() != numBranches) {
            rateCategoryParameter.setDimension(numBranches);
        }

        Parameter randomIndicator = null;
        if (xo.hasChildNamed(RANDOM_INDICATOR)) {
            randomIndicator = (Parameter) xo.getElementFirstChild(RANDOM_INDICATOR);

            if (!randomizeRates) {
                throw new XMLParseException("Cannot provide indicator for randomized rates without randomizeRates=true");
            }

            if (randomIndicator.getDimension() != rateCategoryParameter.getDimension()) {
                throw new XMLParseException("randomIndicator (" + randomIndicator.getDimension()
                        + ") must be the same dimension as the rate parameter (" + rateCategoryParameter.getDimension() + ")");
            }
        }

        Logger.getLogger("dr.evomodel").info("\nUsing an scaled mixture of normals model.");
        Logger.getLogger("dr.evomodel").info("  rates = " + rateCategoryParameter.getDimension());
        Logger.getLogger("dr.evomodel").info("  NB: Make sure you have a prior on "
                + rateCategoryParameter.getId());

        ArbitraryBranchRates.BranchRateTransform transform = parseTransform(xo);

        double scale = xo.getAttribute(RANDOM_SCALE, 1.0);
        if (randomizeRates) {
            for (int i = 0; i < rateCategoryParameter.getDimension(); i++) {
                if (randomIndicator == null || randomIndicator.getParameterValue(i) == 1.0) {
                    double increment = MathUtils.nextGaussian() * scale;
                    double x = transform.randomize(increment);
                    rateCategoryParameter.setValue(i, x);
                }
            }

            if (centerAtOne) {
                double mean = 0.0;
                for (int i = 0; i < rateCategoryParameter.getDimension(); ++i) {
                    mean += rateCategoryParameter.getParameterValue(i);
                }
                mean /= rateCategoryParameter.getDimension();
                
                for (int i = 0; i < rateCategoryParameter.getDimension(); ++i) {
                    rateCategoryParameter.setParameterValue(i,
                            rateCategoryParameter.getParameterValue(i) - mean + 1.0);
                }
            }
        }

        TreeParameterModel.Type includeRoot = xo.getAttribute(INCLUDE_ROOT, false) ?
                TreeParameterModel.Type.WITH_ROOT : TreeParameterModel.Type.WITHOUT_ROOT;

        return new ArbitraryBranchRates(tree, rateCategoryParameter, transform, centerAtOne, includeRoot);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns an arbitrary rate model." +
                "The branch rates are drawn from an arbitrary distribution determine by the prior.";
    }

    static ArbitraryBranchRates.BranchRateTransform parseTransform (XMLObject xo) throws XMLParseException {

        boolean reciprocal = xo.getAttribute(RECIPROCAL, false);
        Logger.getLogger("dr.evomodel").info("  reciprocal = " + reciprocal);

        boolean exp = xo.getAttribute(EXP, false);

        boolean multiplier = xo.getAttribute(MULTIPLIER, false);

        BranchSpecificFixedEffects locationParameter = null;
        if (xo.hasChildNamed(LOCATION)) {
            Object locationObject = xo.getElementFirstChild(LOCATION);
            if ( locationObject instanceof BranchSpecificFixedEffects) {
                locationParameter = (BranchSpecificFixedEffects) locationObject;
            } else {
                locationParameter = new BranchSpecificFixedEffects.None((Parameter) locationObject);
            }
        }

        Parameter scaleParameter = null;
        if (xo.hasChildNamed(SCALE)) {
            scaleParameter = (Parameter) xo.getElementFirstChild(SCALE);
        }

        if (xo.getAttribute(SHRINKAGE, false)) {
            return new ArbitraryBranchRates.BranchRateTransform.LocationShrinkage(
                    ArbitraryBranchRatesParser.ARBITRARY_BRANCH_RATES, locationParameter);
        }

        return make(reciprocal, exp, multiplier, locationParameter, scaleParameter);
    }

    public Class getReturnType() {
        return ArbitraryBranchRates.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class),
            new ElementRule(RATES, Parameter.class, "The rate parameter"),
            AttributeRule.newBooleanRule(RECIPROCAL, true),
            AttributeRule.newBooleanRule(CENTER_AT_ONE, true),
            AttributeRule.newBooleanRule(RANDOMIZE_RATES, true),
            AttributeRule.newBooleanRule(EXP, true),
            AttributeRule.newDoubleRule(RANDOM_SCALE, true),
            new ElementRule(SCALE, Parameter.class, "optional scale parameter", true),
            new ElementRule(LOCATION, Parameter.class, "optional location parameter", true),
            new ElementRule(RANDOM_INDICATOR, new XMLSyntaxRule[] {
                            new ElementRule(Parameter.class),
            }, true),
            AttributeRule.newBooleanRule(SHRINKAGE, true),
    };
}
