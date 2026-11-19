/*
 * CtmcClockRateGradientParser.java
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

package dr.evomodelxml.continuous.hmc;

import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.GradientDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.CtmcClockRateGradient;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.xml.*;

import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;

/**
 * Parses the gradient with respect to the rate of a strict clock, from the cross products of a ctmcRateGradient
 * with the same traitName.
 *
 * @author Marc A. Suchard
 */
public class CtmcClockRateGradientParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "ctmcClockRateGradient";
    private static final String TRAIT_NAME = TreeTraitParserUtilities.TRAIT_NAME;
    private static final String MODE = "mode";

    public String getParserName() {
        return PARSER_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String traitName = xo.getAttribute(TRAIT_NAME, DEFAULT_TRAIT_NAME);
        final TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

        DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        if (!(delegate instanceof GradientDataLikelihoodDelegate)) {
            throw new XMLParseException("Unknown likelihood delegate type");
        }

        BranchRateModel branchRateModel = treeDataLikelihood.getBranchRateModel();
        if (!(branchRateModel instanceof StrictClockBranchRates)) {
            throw new XMLParseException("Only implemented for a strict clock");
        }

        return new CtmcClockRateGradient(traitName, treeDataLikelihood, (GradientDataLikelihoodDelegate) delegate,
                (StrictClockBranchRates) branchRateModel, ApproximateLogCtmcRateGradientParser.parseGradientMode(xo));
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(TRAIT_NAME, true),
            AttributeRule.newStringRule(MODE, true),
            new ElementRule(TreeDataLikelihood.class),
    };

    @Override
    public String getParserDescription() {
        return "Gradient of the log-likelihood with respect to the rate of a strict clock, from the same cross " +
                "products as a ctmcRateGradient";
    }

    @Override
    public Class getReturnType() {
        return CtmcClockRateGradient.class;
    }
}
