/*
 * NodeHeightsStatisticParser.java
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

package dr.evomodelxml.tree;

import dr.evolution.tree.Tree;
import dr.evomodel.tree.NodeHeightsStatistic;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class NodeHeightsStatisticParser extends AbstractXMLObjectParser {

    public static final String NODE_HEIGHTS_STATISTIC = "nodeHeightsStatistic";

    public String getParserName() {
        return NODE_HEIGHTS_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getAttribute(Statistic.NAME, xo.getId());
        Tree tree = (Tree) xo.getChild(Tree.class);

        Parameter groupSizes = (Parameter) xo.getChild(Parameter.class);

        return new NodeHeightsStatistic(name, tree, groupSizes);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A statistic that returns the heights of each internal node in increasing order (or groups them by a group size parameter)";
    }

    public Class getReturnType() {
        return NodeHeightsStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(Statistic.NAME, true),
            new ElementRule(TreeModel.class),
            new ElementRule(Parameter.class, true),
    };

}
