/*
 * CharacterTreePriorParser.java
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

package dr.evomodelxml.substmodel;

import dr.evomodel.substmodel.CharacterTreeModel;
import dr.evomodel.substmodel.CharacterTreePrior;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * Parses a characterTreePrior: the prior on the topology of a sampled character tree.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class CharacterTreePriorParser extends AbstractXMLObjectParser {

    public static final String CHARACTER_TREE_PRIOR = "characterTreePrior";
    public static final String TYPE = "type";
    public static final String C = "c";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final CharacterTreeModel tree = (CharacterTreeModel) xo.getChild(CharacterTreeModel.class);
        final CharacterTreePrior.Type type;
        try {
            type = CharacterTreePrior.Type.parse(xo.getAttribute(TYPE, "uniform"));
        } catch (IllegalArgumentException e) {
            throw new XMLParseException("characterTreePrior: " + e.getMessage());
        }
        if (type == CharacterTreePrior.Type.POLYTOMY && !xo.hasAttribute(C)) {
            throw new XMLParseException("characterTreePrior: type polytomy needs c (c > 1 favours polytomies)");
        }
        if (type != CharacterTreePrior.Type.POLYTOMY && xo.hasAttribute(C)) {
            throw new XMLParseException("characterTreePrior: c applies only to type polytomy");
        }
        final double c = xo.getAttribute(C, 1.0);
        final CharacterTreePrior prior;
        try {
            prior = new CharacterTreePrior(tree, type, c);
        } catch (IllegalArgumentException e) {
            throw new XMLParseException("characterTreePrior: " + e.getMessage());
        }
        if (xo.hasId()) {
            prior.setId(xo.getId());
        }
        final int states = tree.getExternalNodeCount();
        Logger.getLogger("dr.evomodel").info("Character tree prior on '" + tree.getId() + "': " + type.name().toLowerCase()
                + (type == CharacterTreePrior.Type.POLYTOMY ? " (c = " + c + ")" : "") + "; the star tree (F81) has prior "
                + "probability " + Math.exp(prior.logTreeProbability(0)) + " among " + Math.round(Math.exp(prior.getLogTreeCount()))
                + " multifurcating trees on " + states + " states");
        return prior;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(TYPE, true),
            AttributeRule.newDoubleRule(C, true),
            new ElementRule(CharacterTreeModel.class)
    };

    @Override
    public String getParserDescription() {
        return "The prior on the topology of a sampled character tree, corrected for the binary skeleton";
    }

    @Override
    public Class getReturnType() {
        return CharacterTreePrior.class;
    }

    @Override
    public String getParserName() {
        return CHARACTER_TREE_PRIOR;
    }
}
