/*
 * CharacterTreeSPRParser.java
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

import dr.evomodel.operators.CharacterTreeSPR;
import dr.evomodel.substmodel.CharacterTreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 * Parses characterTreeSPR, a root-preserving subtree prune and regraft.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class CharacterTreeSPRParser extends AbstractXMLObjectParser {

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final CharacterTreeModel tree = (CharacterTreeModel) xo.getChild(CharacterTreeModel.class);
        try {
            return new CharacterTreeSPR(tree, xo.getDoubleAttribute(MCMCOperator.WEIGHT));
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
            new ElementRule(CharacterTreeModel.class)
    };

    @Override
    public String getParserDescription() {
        return "A root-preserving subtree prune and regraft on the binary skeleton of a variable character tree";
    }

    @Override
    public Class getReturnType() {
        return CharacterTreeSPR.class;
    }

    @Override
    public String getParserName() {
        return CharacterTreeSPR.CHARACTER_TREE_SPR;
    }
}
