/*
 * CtmcRateGradientParser.java
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

import dr.evomodel.substmodel.LogRateSubstitutionModel;
import dr.evomodel.substmodel.SVSComplexSubstitutionModel;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.GradientDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.AbstractLogAdditiveSubstitutionModelGradient;
import dr.evomodel.treedatalikelihood.discrete.CtmcRateGradient;
import dr.evomodel.treedatalikelihood.discrete.LogCtmcRateGradient;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.xml.*;

import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;

/**
 * Parses the gradient with respect to the rates of a generalSubstitutionModel with a rateIndicator, or with respect
 * to the parameter of a transformedRateSubstitutionModel or logRateSubstitutionModel (the rates are a transform of it).
 *
 * @author Marc A. Suchard
 */
public class CtmcRateGradientParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "ctmcRateGradient";
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

        final LogRateSubstitutionModel transformedModel =
                (LogRateSubstitutionModel) xo.getChild(LogRateSubstitutionModel.class);
        if (transformedModel != null) {
            return new LogCtmcRateGradient(traitName, treeDataLikelihood, (GradientDataLikelihoodDelegate) delegate,
                    transformedModel, ApproximateLogCtmcRateGradientParser.parseGradientMode(xo), false);
        }

        final SVSComplexSubstitutionModel substitutionModel =
                (SVSComplexSubstitutionModel) xo.getChild(SVSComplexSubstitutionModel.class);

        return new CtmcRateGradient(traitName, treeDataLikelihood, (GradientDataLikelihoodDelegate) delegate,
                substitutionModel, ApproximateLogCtmcRateGradientParser.parseGradientMode(xo));
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(TRAIT_NAME, true),
            AttributeRule.newStringRule(MODE, true),
            new ElementRule(TreeDataLikelihood.class),
            new XORRule(
                    new ElementRule(SVSComplexSubstitutionModel.class),
                    new ElementRule(LogRateSubstitutionModel.class)),
    };

    @Override
    public String getParserDescription() {
        return "Gradient of the log-likelihood with respect to the rates of a generalSubstitutionModel " +
                "with a rateIndicator, or to the parameter of a transformedRateSubstitutionModel";
    }

    @Override
    public Class getReturnType() {
        return AbstractLogAdditiveSubstitutionModelGradient.class;
    }
}
