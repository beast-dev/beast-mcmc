/*
 * BayesianBridgeGlobalScaleGradientParser.java
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

import dr.evomodel.branchratemodel.AutoCorrelatedBranchRatesDistribution;
import dr.evomodel.branchratemodel.BayesianBridgeGlobalScaleGradient;
import dr.inferencexml.operators.hmc.HamiltonianMonteCarloOperatorParser;
import dr.xml.*;

/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class BayesianBridgeGlobalScaleGradientParser extends AbstractXMLObjectParser {

    private static final String NAME = "bayesianBridgeGlobalScaleGradient";
    private static final String TOLERANCE = HamiltonianMonteCarloOperatorParser.GRADIENT_CHECK_TOLERANCE;

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        AutoCorrelatedBranchRatesDistribution ratesPrior =
                (AutoCorrelatedBranchRatesDistribution) xo.getChild(AutoCorrelatedBranchRatesDistribution.class);
        try {
            return new BayesianBridgeGlobalScaleGradient(ratesPrior, xo.getAttribute(TOLERANCE, 1E-4));
        } catch (IllegalArgumentException e) {
            throw new XMLParseException(e.getMessage());
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(AutoCorrelatedBranchRatesDistribution.class),
                AttributeRule.newDoubleRule(TOLERANCE, true),
        };
    }

    @Override
    public String getParserDescription() {
        return "Gradient of a Bayesian-bridge auto-correlated rates prior with respect to the bridge's global scale";
    }

    @Override
    public Class getReturnType() {
        return BayesianBridgeGlobalScaleGradient.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}
