/*
 * TreeTimeMeanRateGradientParser.java
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

import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.TreeTimeMeanRateGradient;
import dr.xml.*;

/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class TreeTimeMeanRateGradientParser extends AbstractXMLObjectParser {

    private static final String NAME = "treeTimeMeanRateGradient";
    private static final String TRAIT_NAME = "traitName";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        TreeDataLikelihood likelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
        if (!(likelihood.getDataLikelihoodDelegate() instanceof BeagleDataLikelihoodDelegate)) {
            throw new XMLParseException("Only implemented for sequence (BEAGLE) tree likelihoods");
        }
        BeagleDataLikelihoodDelegate delegate = (BeagleDataLikelihoodDelegate) likelihood.getDataLikelihoodDelegate();
        if (!delegate.isUsePreOrder()) {
            throw new XMLParseException("To use gradients TreeDataLikelihood must have attribute usePreOrder=\"true\"");
        }
        try {
            return new TreeTimeMeanRateGradient(xo.getAttribute(TRAIT_NAME, "Sequence"), likelihood, delegate);
        } catch (IllegalArgumentException e) {
            throw new XMLParseException(e.getMessage());
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(TreeDataLikelihood.class),
                AttributeRule.newStringRule(TRAIT_NAME, true),
        };
    }

    @Override
    public String getParserDescription() {
        return "Gradient of a tree likelihood with respect to the mean rate of a branch-rate model scaled by tree time";
    }

    @Override
    public Class getReturnType() {
        return TreeTimeMeanRateGradient.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}
