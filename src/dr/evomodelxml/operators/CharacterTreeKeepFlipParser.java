/*
 * CharacterTreeKeepFlipParser.java
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

package dr.evomodelxml.operators;

import dr.evomodel.operators.CharacterTreeKeepFlip;
import dr.evomodel.substmodel.CharacterTreeModel;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 * Parses characterTreeKeepFlip, which collapses or restores a node of a sampled character tree.
 *
 * @author Xiang Ji
 */
public class CharacterTreeKeepFlipParser extends AbstractXMLObjectParser {

    public static final String RATES = "rates";
    public static final String SLAB = "slab";
    public static final String NEAR_ZERO_WEIGHT = "nearZeroWeight";
    public static final String NEAR_ZERO_SCALE = "nearZeroScale";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final CharacterTreeModel tree = (CharacterTreeModel) xo.getChild(CharacterTreeModel.class);
        final Parameter rates = (Parameter) xo.getElementFirstChild(RATES);
        final ParametricDistributionModel slab = (ParametricDistributionModel) xo.getElementFirstChild(SLAB);
        try {
            return new CharacterTreeKeepFlip(tree, rates, slab, xo.getAttribute(NEAR_ZERO_WEIGHT, 0.5),
                    xo.getAttribute(NEAR_ZERO_SCALE, 0.1), xo.getDoubleAttribute(MCMCOperator.WEIGHT));
        } catch (IllegalArgumentException e) {
            throw new XMLParseException(getParserName() + ": " + e.getMessage());
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newDoubleRule(NEAR_ZERO_WEIGHT, true),
            AttributeRule.newDoubleRule(NEAR_ZERO_SCALE, true),
            new ElementRule(CharacterTreeModel.class),
            new ElementRule(RATES, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(SLAB, new XMLSyntaxRule[]{new ElementRule(ParametricDistributionModel.class)})
    };

    @Override
    public String getParserDescription() {
        return "Collapses or restores one non-root node of a sampled character tree, drawing a new rate";
    }

    @Override
    public Class getReturnType() {
        return CharacterTreeKeepFlip.class;
    }

    @Override
    public String getParserName() {
        return CharacterTreeKeepFlip.CHARACTER_TREE_KEEP_FLIP;
    }
}
