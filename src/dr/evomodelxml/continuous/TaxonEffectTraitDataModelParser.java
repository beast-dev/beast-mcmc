/*
 * TaxonEffectTraitDataModelParser.java
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

package dr.evomodelxml.continuous;

import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.continuous.TaxonEffectTraitDataModel;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.*;

import static dr.evomodelxml.continuous.ContinuousTraitDataModelParser.NUM_TRAITS;

public class TaxonEffectTraitDataModelParser extends AbstractXMLObjectParser {

    private static final String CONTINUOUS_TRAITS = "taxonEffectTraitDataModel";
    private static final String TAXON_EFFECTS = "effects";
    private static final String SIGN = "sign";
    private static final String SET_NAMES = "setEffectParameterNames";
    private static final String CHECK_NAMES = "checkEffectParameterNames";
    private static final String RANDOMIZE = "randomizeEffects";
    private static final String RANDOMIZATION_ST_DEV = "randomizationStDev";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        PrecisionType precisionType = PrecisionType.SCALAR;

        Tree treeModel = (Tree) xo.getChild(Tree.class);
        TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();
        TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
                utilities.parseTraitsFromTaxonAttributes(xo, treeModel, true);
        CompoundParameter traitParameter = returnValue.traitParameter;

        int dimAll = traitParameter.getParameter(0).getDimension();
        int numTraits = xo.getAttribute(NUM_TRAITS, 1);
        int dim = dimAll / numTraits;

        boolean[] missingIndicators = returnValue.getMissingIndicators();
        String traitName = returnValue.traitName;
        boolean useMissingIndices = returnValue.useMissingIndices;

        if (xo.hasChildNamed(TreeTraitParserUtilities.JITTER)) {
            utilities.jitter(xo, dim, missingIndicators);
        }

        final TaxonEffectTraitDataModel.EffectMap map;
        final Parameter effects;

        if (xo.hasChildNamed(TAXON_EFFECTS)) {
            effects = (Parameter) xo.getElementFirstChild(TAXON_EFFECTS);
            Parameter sign = null;
            if (xo.hasChildNamed(SIGN)) {
                sign = (Parameter) xo.getElementFirstChild(SIGN);
            }
            map = new TaxonEffectTraitDataModel.EffectMap(treeModel, effects, sign, dim);
        } else {
            TaxonEffectTraitDataModel original = (TaxonEffectTraitDataModel)
                    xo.getChild(TaxonEffectTraitDataModel.class);
            map = new TaxonEffectTraitDataModel.EffectMap(treeModel, original.getMap());
            effects = map.getEffects();
        }

        if (effects.getDimension() != traitParameter.getDimension()) {
            throw new XMLParseException("Invalid effects dimension");
        }

        if (xo.getAttribute(RANDOMIZE, false)) {
            double stdDev = xo.getAttribute(RANDOMIZATION_ST_DEV, 1.0);
            for (int i = 0; i < effects.getDimension(); ++i) {
                effects.setParameterValue(i, effects.getParameterValue(i) +
                        MathUtils.nextGaussian() * stdDev);
            }
        }

        if (treeModel.getExternalNodeCount() != effects.getDimension() * dim) {
            throw new XMLParseException("Invalid effect dimension");
        }

        return new TaxonEffectTraitDataModel(traitName,
                treeModel,
                traitParameter,
                map,
                missingIndicators, useMissingIndices,
                dim, precisionType);
    }

    public static final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(Tree.class),
            new ElementRule(TreeTraitParserUtilities.TRAIT_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new XORRule(
                    new ElementRule(TAXON_EFFECTS, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }),
                    new ElementRule(TaxonEffectTraitDataModel.class)),
            new ElementRule(SIGN, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class),
            }, true),
            AttributeRule.newIntegerRule(NUM_TRAITS, true),
            AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME, true),
            AttributeRule.newBooleanRule(SET_NAMES, true),
            AttributeRule.newBooleanRule(CHECK_NAMES, true),
            AttributeRule.newBooleanRule(RANDOMIZE, true),
            AttributeRule.newDoubleRule(RANDOMIZATION_ST_DEV, true),
    };

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return "parses continuous traits with taxon effects from a tree";
    }

    @Override
    public Class getReturnType() {
        return TaxonEffectTraitDataModel.class;
    }

    @Override
    public String getParserName() {
        return CONTINUOUS_TRAITS;
    }
}
