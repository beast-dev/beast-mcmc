/*
 * AutoCorrelatedGradientWrtIncrementsParser.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evomodel.branchratemodel.AutoCorrelatedBranchRatesDistribution;
import dr.evomodel.branchratemodel.AutoCorrelatedGradientWrtIncrements;
import dr.evomodel.branchratemodel.BranchRateGradientWrtIncrements;
import dr.evomodel.treedatalikelihood.continuous.BranchRateGradient;
import dr.evomodel.treedatalikelihood.discrete.BranchRateGradientForDiscreteTrait;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.xml.*;

public class BranchRateGradientWrtIncrementsParser extends AbstractXMLObjectParser {

    private static final String GRADIENT = "branchRateGradientWrtIncrements";

    public String getParserName() {
        return GRADIENT;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        AutoCorrelatedGradientWrtIncrements priorProvider = (AutoCorrelatedGradientWrtIncrements)
                xo.getChild(AutoCorrelatedGradientWrtIncrements.class);

        GradientWrtParameterProvider rateProvider = (GradientWrtParameterProvider)
                xo.getChild(GradientWrtParameterProvider.class);

        if (!(rateProvider instanceof BranchRateGradient)) {
            throw new XMLParseException("Must provide a branch rate gradient");
        }

        return new BranchRateGradientWrtIncrements(rateProvider, priorProvider);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns an arbitrary rate model." +
                "The branch rates are drawn from an arbitrary distribution determine by the prior.";
    }

    public Class getReturnType() {
        return AutoCorrelatedGradientWrtIncrements.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(AutoCorrelatedGradientWrtIncrements.class),
            new XORRule(
                new ElementRule(BranchRateGradient.class),
                new ElementRule(BranchRateGradientForDiscreteTrait.class)
            ),
    };
}
