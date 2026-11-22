/*
 * TaxonEffectTraitDataModelParser.java
 *
 * Copyright © 2002-2025 the BEAST Development Team
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
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

public class TreeTraitAttributeParser extends AbstractXMLObjectParser {

    private static final String TRAIT_ATTRIBUTE_PARSER = "treeTraitAttribute";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Tree treeModel = (Tree) xo.getChild(Tree.class);
        TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();
        TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
                utilities.parseTraitsFromTaxonAttributes(xo, treeModel, true);
        CompoundParameter traitParameter = returnValue.traitParameter;

        int dimAll = traitParameter.getParameter(0).getDimension();
        int numTraits = xo.getAttribute(ContinuousTraitDataModelParser.NUM_TRAITS, 1);
        int dim = dimAll / numTraits;

        boolean[] missingIndicators = returnValue.getMissingIndicators();

        if (xo.hasChildNamed(TreeTraitParserUtilities.JITTER)) {
            utilities.jitter(xo, dim, missingIndicators);
        }

        return null;
    }

    public static final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(Tree.class),
            new ElementRule(TreeTraitParserUtilities.TRAIT_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            AttributeRule.newIntegerRule(ContinuousTraitDataModelParser.NUM_TRAITS, true),
            AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME, true),
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
        return null;
    }

    @Override
    public String getParserName() {
        return TRAIT_ATTRIBUTE_PARSER;
    }
}
