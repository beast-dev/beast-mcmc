/*
 * FastTreeIntervalsParser.java
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

package dr.evomodelxml.coalescent.intervals;

import dr.evomodel.coalescent.intervals.FastTreeIntervals;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;


public class FastTreeIntervalsParser extends AbstractXMLObjectParser {

    public static final String TREE_INTERVALS = "fastTreeIntervals";
    // the name used before BigFastTreeIntervals was renamed to FastTreeIntervals
    public static final String OLD_TREE_INTERVALS = "bigFastTreeIntervals";
    public static final String TREE = "tree";

    public String getParserName() {
        return TREE_INTERVALS;
    }

    public String[] getParserNames() {
        return new String[]{TREE_INTERVALS, OLD_TREE_INTERVALS};
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        return new FastTreeIntervals( tree);

    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Parser for fastTreeIntervals.";
    }

    public Class getReturnType() {
        return FastTreeIntervals.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class)
    };

}