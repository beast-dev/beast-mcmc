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

import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.AutoCorrelatedBranchRatesDistribution;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DifferentiableBranchRatesFullMethods;
import dr.inference.distribution.ParametricMultivariateDistributionModel;
import dr.xml.*;

/**
 */
public class AutoCorrelatedBranchRatesDistributionParser extends AbstractXMLObjectParser {

    private static final String AUTO_CORRELATED_RATES = "autoCorrelatedRatesPrior";
    private static final String SCALING = "scaling";
    private static final String LOG = "log";

    public String getParserName() {
        return AUTO_CORRELATED_RATES;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        DifferentiableBranchRatesFullMethods branchRates = (DifferentiableBranchRatesFullMethods) xo.getChild(BranchRateModel.class);
        if(!(branchRates instanceof DifferentiableBranchRatesFullMethods)){
            throw new XMLParseException("Branch rate model must be implement DifferentialBranchRateMethodsFul");
        }

        ParametricMultivariateDistributionModel distribution = (ParametricMultivariateDistributionModel)
                           xo.getChild(ParametricMultivariateDistributionModel.class);

        AutoCorrelatedBranchRatesDistribution.BranchVarianceScaling scaling = parseScaling(xo);

        boolean log = xo.getAttribute(LOG, false);

        return new AutoCorrelatedBranchRatesDistribution(xo.getId(), branchRates, distribution,
                scaling, log);
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
            new ElementRule(BranchRateModel.class),
            new ElementRule(ParametricMultivariateDistributionModel.class),
            AttributeRule.newStringRule(SCALING, true),
            AttributeRule.newBooleanRule(LOG, true),
    };
}
