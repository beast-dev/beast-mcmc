/*
 * BranchRateGradientWrtIncrementsParser.java
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

import dr.evomodel.branchratemodel.AutoCorrelatedGradientWrtIncrements;
import dr.evomodel.branchratemodel.BranchRateGradientWrtIncrements;
import dr.evomodel.treedatalikelihood.continuous.BranchRateGradient;
import dr.evomodel.treedatalikelihood.continuous.BranchSpecificOptimaGradient;
import dr.evomodel.treedatalikelihood.discrete.BranchRateGradientForDiscreteTrait;
import dr.inference.hmc.CompoundGradient;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.JointBranchRateGradient;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

public class BranchRateGradientWrtIncrementsParser extends AbstractXMLObjectParser {

    private static final String GRADIENT = "branchRateGradientWrtIncrements";

    public String getParserName() {
        return GRADIENT;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {


        List<AutoCorrelatedGradientWrtIncrements> priorProvider = new ArrayList<AutoCorrelatedGradientWrtIncrements>();

        AutoCorrelatedGradientWrtIncrements onePriorProvider = (AutoCorrelatedGradientWrtIncrements)
                xo.getChild(AutoCorrelatedGradientWrtIncrements.class);
        if (onePriorProvider != null) {
            priorProvider.add(onePriorProvider);
        } else {
            CompoundGradient cxo = (CompoundGradient) xo.getChild(CompoundGradient.class);
            for (GradientWrtParameterProvider gradientChild : cxo.getDerivativeList()) {
                if (gradientChild instanceof AutoCorrelatedGradientWrtIncrements) {
                    priorProvider.add((AutoCorrelatedGradientWrtIncrements) gradientChild);
                } else {
                    throw new XMLParseException("Compound gradient must only contain AutoCorrelatedGradientWrtIncrements and not " + (gradientChild.getClass()));
                }
            }
        }

        //        BranchSpecificOptimaGradient rateProvider = (BranchSpecificOptimaGradient) xo.getChild(BranchSpecificOptimaGradient.class);

        GradientWrtParameterProvider rateProvider = null;

        for (GradientWrtParameterProvider gradientProvider : xo.getAllChildren(GradientWrtParameterProvider.class)) {
            if (!(gradientProvider instanceof AutoCorrelatedGradientWrtIncrements) && !(gradientProvider instanceof CompoundGradient)) {
                rateProvider = gradientProvider;
            }
        }

        if (!(rateProvider instanceof JointBranchRateGradient) && !(rateProvider instanceof BranchSpecificOptimaGradient)) {
            if (!(rateProvider instanceof BranchRateGradient) &&
                    !(rateProvider instanceof BranchRateGradientForDiscreteTrait)) {
                throw new XMLParseException("Must provide a branch rate gradient");
            }
        }
        return new BranchRateGradientWrtIncrements(rateProvider, priorProvider);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns the gradient of an arbitrary branch rate model " +
                "w.r.t. the increments in an auto-correlated rates prior.";
    }

    public Class getReturnType() {
        return AutoCorrelatedGradientWrtIncrements.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new XORRule(
                    new ElementRule(AutoCorrelatedGradientWrtIncrements.class),
                    new ElementRule(CompoundGradient.class)
            ),
            new XORRule(
                    new XORRule(
                            new ElementRule(BranchRateGradient.class),
                            new ElementRule(BranchSpecificOptimaGradient.class)
                    ),
                    new XORRule(
                            new ElementRule(BranchRateGradientForDiscreteTrait.class),
                            new ElementRule(JointBranchRateGradient.class)
                    )
            ),
    };
}
