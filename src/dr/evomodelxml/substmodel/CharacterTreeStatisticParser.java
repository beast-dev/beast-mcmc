/*
 * CharacterTreeStatisticParser.java
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
import dr.evomodel.substmodel.CharacterTreeStatistic;
import dr.evomodel.substmodel.HierarchicalSubstitutionModel;
import dr.xml.*;

import java.util.List;
import java.util.logging.Logger;

/**
 * Parses a characterTreeStatistic; at S <= 5 it prints the list of trees that the model column indexes.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class CharacterTreeStatisticParser extends AbstractXMLObjectParser {

    public static final String CHARACTER_TREE_STATISTIC = "characterTreeStatistic";
    public static final String CLADES = "clades";
    public static final String EXCHANGEABILITIES = "exchangeabilities";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final CharacterTreeModel tree = (CharacterTreeModel) xo.getChild(CharacterTreeModel.class);
        final HierarchicalSubstitutionModel model = (HierarchicalSubstitutionModel) xo.getChild(HierarchicalSubstitutionModel.class);
        final int states = tree.getExternalNodeCount();
        final boolean clades = xo.getAttribute(CLADES, states <= CharacterTreeStatistic.MAX_STATES_FOR_CLADES);
        final boolean exchangeabilities = xo.getAttribute(EXCHANGEABILITIES,
                model != null && states <= CharacterTreeStatistic.MAX_STATES_FOR_EXCHANGEABILITIES);
        final String id = xo.hasId() ? xo.getId() : tree.getId();
        final CharacterTreeStatistic statistic;
        try {
            statistic = new CharacterTreeStatistic(id, tree, model, clades, exchangeabilities);
        } catch (IllegalArgumentException e) {
            throw new XMLParseException("characterTreeStatistic '" + id + "': " + e.getMessage());
        }
        final List<String> trees = statistic.getTreeList();
        if (trees != null) {
            final StringBuilder sb = new StringBuilder("Character tree statistic '" + id + "': column " + id
                    + ".model indexes these " + trees.size() + " trees:");
            for (int i = 0; i < trees.size(); i++) {
                sb.append("\n  ").append(i).append(' ').append(trees.get(i));
            }
            Logger.getLogger("dr.evomodel").info(sb.toString());
        }
        return statistic;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(CLADES, true),
            AttributeRule.newBooleanRule(EXCHANGEABILITIES, true),
            new ElementRule(CharacterTreeModel.class),
            new ElementRule(HierarchicalSubstitutionModel.class, true)
    };

    @Override
    public String getParserDescription() {
        return "The number of clades, the star indicator, the tree index, clade indicators and exchangeabilities of a "
                + "sampled character tree";
    }

    @Override
    public Class getReturnType() {
        return CharacterTreeStatistic.class;
    }

    @Override
    public String getParserName() {
        return CHARACTER_TREE_STATISTIC;
    }
}
