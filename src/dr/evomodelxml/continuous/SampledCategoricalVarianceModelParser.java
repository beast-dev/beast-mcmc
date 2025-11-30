/*
 * ContinuousTraitDataModelParser.java
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
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitDataModel;
import dr.evomodelxml.continuous.dr.evomodel.continuous.SampledCategoricalVarianceModel;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.*;

public class SampledCategoricalVarianceModelParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "sampledCategoricalVarianceModel";
    private static final String TRAITS = "traitParameter";
    private static final String VARIANCES = "variances";
    private static final String ASSIGNMENTS = "assignments";
    private static final String RANDOMIZE_TRAITS = "randomizeTraits";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Tree treeModel = (Tree) xo.getChild(Tree.class);
        TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();

        TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
                utilities.parseTraitsFromTaxonAttributes(xo, treeModel, true);
        CompoundParameter oldMeans = returnValue.traitParameter;
        Parameter means = new Parameter.Default(oldMeans.getId() + ".mean", oldMeans.getParameterValues());

        Parameter traits = (Parameter) xo.getElementFirstChild(TRAITS);
        Parameter variances = (Parameter) xo.getElementFirstChild(VARIANCES);
        Parameter assignments = (Parameter) xo.getElementFirstChild(ASSIGNMENTS);

        if (xo.getAttribute(RANDOMIZE_TRAITS, false)) {
            double scale = 1.0;
            for (int i = 0; i < traits.getDimension(); ++i) {
                traits.setParameterValueQuietly(i,
                        traits.getParameterValue(i) + scale * MathUtils.nextGaussian());
            }
            traits.fireParameterChangedEvent();
        }

        return new SampledCategoricalVarianceModel(xo.getId(), traits, means, variances, assignments);
    }

    private static final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(Tree.class),
            new ElementRule(TreeTraitParserUtilities.TRAIT_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME, true),
            new ElementRule(TRAITS, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(VARIANCES, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(ASSIGNMENTS, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            AttributeRule.newBooleanRule(RANDOMIZE_TRAITS, true),
    };

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return "parses continuous traits from a tree";
    }

    @Override
    public Class getReturnType() {
        return ContinuousTraitDataModel.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }
}
