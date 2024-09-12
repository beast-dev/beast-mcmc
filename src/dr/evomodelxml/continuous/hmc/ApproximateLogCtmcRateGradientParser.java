/*
 * ApproximateLogCtmcRateGradientParser.java
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

import dr.evomodel.substmodel.GlmSubstitutionModel;
import dr.evomodel.substmodel.LogRateSubstitutionModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.AbstractGlmSubstitutionModelGradient;
import dr.evomodel.treedatalikelihood.discrete.LogCtmcRateGradient;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.xml.*;

import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;

/**
 * @author Filippo Monti
 * @author Marc A. Suchard
 */

public class ApproximateLogCtmcRateGradientParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "approximateLogCtmcRateGradient";
    private static final String TRAIT_NAME = TreeTraitParserUtilities.TRAIT_NAME;

    public String getParserName(){ return PARSER_NAME; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String traitName = xo.getAttribute(TRAIT_NAME, DEFAULT_TRAIT_NAME);
        final TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

        DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        if (!(delegate instanceof BeagleDataLikelihoodDelegate)) {
            throw new XMLParseException("Unknown likelihood delegate type");
        }

        for (Object child : xo.getChildren()) {
            if (child instanceof GlmSubstitutionModel) {
                final GlmSubstitutionModel substitutionModel = (GlmSubstitutionModel) child;
                return new LogCtmcRateGradient(traitName, treeDataLikelihood,
                        (BeagleDataLikelihoodDelegate) delegate, substitutionModel);
            } else if (child instanceof LogRateSubstitutionModel) {
                final LogRateSubstitutionModel substitutionModel = (LogRateSubstitutionModel) child;
                return new LogCtmcRateGradient(traitName, treeDataLikelihood,
                        (BeagleDataLikelihoodDelegate) delegate, substitutionModel);
            }
        }
        throw new XMLParseException("No valid substitution model found"); // TODO ugly (!?)
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(TRAIT_NAME, true),
            new ElementRule(TreeDataLikelihood.class),
            new XORRule(
                    new ElementRule(GlmSubstitutionModel.class),
                    new ElementRule(LogRateSubstitutionModel.class)),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return AbstractGlmSubstitutionModelGradient.class;
    }
}

