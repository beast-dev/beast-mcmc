/*
 * ArbitraryBranchRatesParser.java
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

import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchSpecificFixedEffects;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.logging.Logger;

import static dr.evomodel.branchratemodel.ArbitraryBranchRates.make;

/**
 */
public class ArbitraryBranchRatesParser extends AbstractXMLObjectParser {

    public static final String ARBITRARY_BRANCH_RATES = "arbitraryBranchRates";
    private static final String RATES = "rates";
    static final String RECIPROCAL = "reciprocal";
    static final String EXP = "exp";
    private static final String MULTIPLIER = "multiplier";
    private static final String CENTER_AT_ONE = "centerAtOne";
    private static final String RANDOMIZE_RATES = "randomizeRates";
    private static final String RANDOM_SCALE = "randomScale";
    static final String LOCATION = "location";
    static final String SCALE = "scale";

    public String getParserName() {
        return ARBITRARY_BRANCH_RATES;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

        XMLObject cxo = xo.getChild(RATES);

        Parameter rateCategoryParameter = (Parameter) cxo.getChild(Parameter.class);

        boolean centerAtOne = xo.getAttribute(CENTER_AT_ONE, true);

        boolean randomizeRates = xo.getAttribute(RANDOMIZE_RATES, false);

        if(centerAtOne && randomizeRates) {
            throw new XMLParseException("Cannot centerAtOne and randomize the starting rates");
        }

        final int numBranches = tree.getNodeCount() - 1;
        if (rateCategoryParameter.getDimension() > 1 && (rateCategoryParameter.getDimension() != numBranches)) {
            throw new XMLParseException("Incorrect number of rate parameters");
        }

        if (rateCategoryParameter.getDimension() != numBranches) {
            rateCategoryParameter.setDimension(numBranches);
        }

        Logger.getLogger("dr.evomodel").info("\nUsing an scaled mixture of normals model.");
        Logger.getLogger("dr.evomodel").info("  rates = " + rateCategoryParameter.getDimension());
        Logger.getLogger("dr.evomodel").info("  NB: Make sure you have a prior on "
                + rateCategoryParameter.getId());


        ArbitraryBranchRates.BranchRateTransform transform = parseTransform(xo);

        double scale = xo.getAttribute(RANDOM_SCALE, 1.0);
        if (randomizeRates) {
            for (int i = 0; i < rateCategoryParameter.getDimension(); i++) {
                double increment = MathUtils.nextGaussian() * scale;
                double x = transform.randomize(increment);
                rateCategoryParameter.setValue(i, x);
            }
        }

        return new ArbitraryBranchRates(tree, rateCategoryParameter, transform, centerAtOne);
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
    };
}
