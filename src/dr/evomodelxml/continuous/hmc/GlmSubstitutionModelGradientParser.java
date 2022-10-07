/*
 * GlmSubstitutionModelGradientParser.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.continuous.hmc;

import dr.evomodel.substmodel.OldGLMSubstitutionModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.DesignMatrixSubstitutionModelGradient;
import dr.evomodel.treedatalikelihood.discrete.AbstractGlmSubstitutionModelGradient;
import dr.evomodel.treedatalikelihood.discrete.FixedEffectSubstitutionModelGradient;
import dr.evomodel.treedatalikelihood.discrete.RandomEffectsSubstitutionModelGradient;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.DesignMatrix;
import dr.inference.model.MaskedParameter;
import dr.xml.*;

import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;

/**
 * @author Marc A. Suchard
 */

public class GlmSubstitutionModelGradientParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "glmSubstitutionModelGradient";
    private static final String TRAIT_NAME = TreeTraitParserUtilities.TRAIT_NAME;
    private static final String EFFECTS = "effects";

    public String getParserName(){ return PARSER_NAME; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String traitName = xo.getAttribute(TRAIT_NAME, DEFAULT_TRAIT_NAME);
        final TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
        OldGLMSubstitutionModel substitutionModel = (OldGLMSubstitutionModel) xo.getChild(OldGLMSubstitutionModel.class);

        DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        if (!(delegate instanceof BeagleDataLikelihoodDelegate)) {
            throw new XMLParseException("Unknown likelihood delegate type");
        }

        String effectsString = xo.getAttribute(EFFECTS, "fixed");
        if (effectsString.equalsIgnoreCase("fixed")) {
            if (substitutionModel.getGeneralizedLinearModel().getNumberOfFixedEffects() < 1) {
                throw new XMLParseException("No fixed effects in '" + substitutionModel.getId() + "'");
            }
            return new FixedEffectSubstitutionModelGradient(traitName, treeDataLikelihood,
                    (BeagleDataLikelihoodDelegate) delegate, substitutionModel) {
            };
        } else if (effectsString.equalsIgnoreCase("random")) {
            if (substitutionModel.getGeneralizedLinearModel().getNumberOfRandomEffects() < 1) {
                throw new XMLParseException("No random effects in '" + substitutionModel.getId() + "'");
            }
            return new RandomEffectsSubstitutionModelGradient(traitName, treeDataLikelihood,
                    (BeagleDataLikelihoodDelegate) delegate, substitutionModel);
        } else if (effectsString.equalsIgnoreCase("design")) {
            MaskedParameter parameter = (MaskedParameter) xo.getChild(MaskedParameter.class);
            DesignMatrix matrix = (DesignMatrix) xo.getChild(DesignMatrix.class);
            return new DesignMatrixSubstitutionModelGradient(traitName, treeDataLikelihood,
                    (BeagleDataLikelihoodDelegate) delegate, substitutionModel, matrix, parameter);
        } else {
            throw new XMLParseException("Unknown effects type '" + effectsString + "'");
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(TRAIT_NAME, true),
            new ElementRule(TreeDataLikelihood.class),
            new ElementRule(OldGLMSubstitutionModel.class),
            AttributeRule.newStringRule(EFFECTS, true),
            new ElementRule(MaskedParameter.class, true),
            new ElementRule(DesignMatrix.class, true),
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
