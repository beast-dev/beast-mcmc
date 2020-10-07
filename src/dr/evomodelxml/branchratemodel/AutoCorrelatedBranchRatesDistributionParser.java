/*
 * AutoCorrelatedBranchRatesDistributionParser.java
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

import dr.evomodel.branchratemodel.*;
import dr.inference.distribution.ParametricMultivariateDistributionModel;
import dr.xml.*;

/**
 */
public class AutoCorrelatedBranchRatesDistributionParser extends AbstractXMLObjectParser {

    private static final String AUTO_CORRELATED_RATES = "autoCorrelatedRatesPrior";
    private static final String SCALING = "scaling";
    private static final String LOG = "log";
    private static final String OPERATE_ON_INCREMENTS = "operateOnIncrements";

    public String getParserName() {
        return AUTO_CORRELATED_RATES;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        DifferentiableBranchRates branchRates = (DifferentiableBranchRates) xo.getChild(BranchRateModel.class);

        ParametricMultivariateDistributionModel distribution = (ParametricMultivariateDistributionModel)
                           xo.getChild(ParametricMultivariateDistributionModel.class);

        AutoCorrelatedBranchRatesDistribution.BranchVarianceScaling scaling = parseScaling(xo);

        boolean log = xo.getAttribute(LOG, false);

        boolean operateOnIncrements = xo.getAttribute(OPERATE_ON_INCREMENTS, false);

        // TODO Change parser to accept Tree and then pass to ACBRD
        return new AutoCorrelatedBranchRatesDistribution(xo.getId(), branchRates, distribution,
                scaling, log, operateOnIncrements);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns an arbitrary rate model." +
                "The branch rates are drawn from an arbitrary distribution determine by the prior.";
    }

    private static AutoCorrelatedBranchRatesDistribution.BranchVarianceScaling parseScaling(XMLObject xo)
            throws XMLParseException {

        String defaultName = AutoCorrelatedBranchRatesDistribution.BranchVarianceScaling.NONE.getName();
        String name = xo.getAttribute(SCALING, defaultName);

        AutoCorrelatedBranchRatesDistribution.BranchVarianceScaling scaling =
                AutoCorrelatedBranchRatesDistribution.BranchVarianceScaling.parse(name);

        if (scaling == null) {
            throw new XMLParseException("Unknown scaling type");
        }

        return scaling;
    }

    public Class getReturnType() {
        return ArbitraryBranchRates.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(DifferentiableBranchRates.class),
            new ElementRule(ParametricMultivariateDistributionModel.class),
            AttributeRule.newStringRule(SCALING, true),
            AttributeRule.newBooleanRule(LOG, true),
            AttributeRule.newBooleanRule(OPERATE_ON_INCREMENTS, true),
    };
}
