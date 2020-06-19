/*
 * TreeHeightStatisticParser.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.tree;

import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeHeightStatistic;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 * Joseph Heled
 */
public class TreeHeightStatisticParser extends AbstractXMLObjectParser {

    public static final String TREE_HEIGHT_STATISTIC = "treeHeightStatistic";

        public String getParserName() {
            return TREE_HEIGHT_STATISTIC;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            final String name = xo.getAttribute(Statistic.NAME, xo.getId());
            final Tree tree = (Tree) xo.getChild(Tree.class);

            return new TreeHeightStatistic(name, tree);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A statistic that returns the height of the tree";
        }

        public Class getReturnType() {
            return TreeHeightStatistic.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule(Statistic.NAME, true),
                new ElementRule(Tree.class),
        };
}