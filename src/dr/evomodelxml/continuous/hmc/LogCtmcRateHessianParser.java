/*
 * LogCtmcRateHessianParser.java
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
import dr.evomodel.substmodel.LogAdditiveCtmcRateProvider;
import dr.evomodel.substmodel.LogRateSubstitutionModel;
import dr.evomodel.substmodel.StronglyLumpableCtmcRates;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.GradientDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.AbstractLogAdditiveSubstitutionModelGradient;
import dr.evomodel.treedatalikelihood.discrete.LogCtmcRateHessian;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;
import dr.xml.XORRule;

import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;

public class LogCtmcRateHessianParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "logCtmcRateHessian";
    private static final String APPROXIMATE_PARSER_NAME = "approximateLogCtmcRateHessian";
    private static final String EXACT_PARSER_NAME = "exactLogCtmcRateHessian";
    private static final String TRAIT_NAME = TreeTraitParserUtilities.TRAIT_NAME;
    private static final String FORCE_ALL_REAL = "forceAllReal";
    private static final String MODE = "mode";

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    @Override
    public String[] getParserNames() {
        return new String[] { PARSER_NAME, APPROXIMATE_PARSER_NAME, EXACT_PARSER_NAME };
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        String traitName = xo.getAttribute(TRAIT_NAME, DEFAULT_TRAIT_NAME);
        final boolean forceAllReal = xo.getAttribute(FORCE_ALL_REAL, false);
        final TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

        DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        if (!(delegate instanceof GradientDataLikelihoodDelegate)) {
            throw new XMLParseException("Unknown likelihood delegate type");
        }

        for (Object child : xo.getChildren()) {
            if (child instanceof GlmSubstitutionModel) {
                if (!(delegate instanceof BeagleDataLikelihoodDelegate)) {
                    throw new XMLParseException("GLM CTMC Hessian requires a BeagleDataLikelihoodDelegate");
                }
                return new LogCtmcRateHessian(traitName, treeDataLikelihood,
                        (BeagleDataLikelihoodDelegate) delegate, (GlmSubstitutionModel) child);
            } else if (child instanceof LogRateSubstitutionModel) {
                final LogRateSubstitutionModel substitutionModel = (LogRateSubstitutionModel) child;
                LogAdditiveCtmcRateProvider rates = substitutionModel.getRateProvider();
                if (rates instanceof StronglyLumpableCtmcRates) {
                    throw new XMLParseException("Log CTMC rate Hessian is not implemented for strongly lumpable rates");
                }

                AbstractLogAdditiveSubstitutionModelGradient.ApproximationMode mode = parseMode(xo);
                if (xo.getName().equalsIgnoreCase(EXACT_PARSER_NAME)) {
                    mode = AbstractLogAdditiveSubstitutionModelGradient.ApproximationMode.EXACT_SPECTRAL;
                }

                return new LogCtmcRateHessian(traitName, treeDataLikelihood,
                        (GradientDataLikelihoodDelegate) delegate, substitutionModel, mode, forceAllReal);
            }
        }

        throw new XMLParseException("No valid substitution model found");
    }

    private AbstractLogAdditiveSubstitutionModelGradient.ApproximationMode parseMode(XMLObject xo)
            throws XMLParseException {
        String modeString = xo.getAttribute(MODE, "firstOrder");
        if (modeString.compareToIgnoreCase("exact") == 0) {
            return AbstractLogAdditiveSubstitutionModelGradient.ApproximationMode.EXACT_SPECTRAL;
        }
        try {
            return AbstractLogAdditiveSubstitutionModelGradient.ApproximationMode.factory(modeString);
        } catch (Exception e) {
            throw new XMLParseException(e.getMessage());
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(TRAIT_NAME, true),
            AttributeRule.newBooleanRule(FORCE_ALL_REAL, true),
            AttributeRule.newStringRule(MODE, true),
            new ElementRule(TreeDataLikelihood.class),
            new XORRule(
                    new XMLSyntaxRule[]{
                            new ElementRule(GlmSubstitutionModel.class),
                            new ElementRule(LogRateSubstitutionModel.class),
                    }),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return LogCtmcRateHessian.class;
    }
}
