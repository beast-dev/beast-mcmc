/*
 * LineageCountThroughTimeParser.java
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

package dr.evoxml;

import dr.evolution.tree.LineageCountThroughTime;
import dr.stats.Variate;
import dr.xml.*;

/**
 * @author Andrew Rambaut
 * @version $Id: DateParser.java,v 1.2 2005/05/24 20:25:59 rambaut Exp $
 */
public class LineageCountThroughTimeParser extends AbstractXMLObjectParser {

    public static final String TREE_FILE = "treeFile";
    public static final String MIN_TIME = "minTime";
    public static final String MAX_TIME = "maxTime";
    public static final String BIN_COUNT = "binCount";
    public static final String SKIP = "skip";

    public String getParserName() {
        return "lineageCountThroughTime";
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        if (xo.getChildCount() > 0) {
            throw new XMLParseException("No child elements allowed in location element.");
        }

        int binCount = xo.getAttribute(BIN_COUNT, 100);
        int skip = xo.getAttribute(SKIP, 0);
        double minTime = xo.getAttribute(MIN_TIME, 0.0);
        double maxTime = xo.getAttribute(MAX_TIME, 1.0);

        String treeFile = xo.getAttribute(TREE_FILE, "");

        try {
            Variate[] ltt = LineageCountThroughTime.getLTT(treeFile, minTime, maxTime, binCount, skip);

            StringBuilder builder = new StringBuilder();
            builder.append("Time\tLineage_Count\n");
            for (int i = 0; i < ltt[0].getCount(); i++) {
                builder.append(ltt[0].get(i));
                builder.append("\t");
                builder.append(ltt[1].get(i));
                builder.append("\n");
            }
            String lttString = builder.toString();

            System.out.println();
            System.out.println(lttString);
            return lttString;

        } catch (Exception e) {

            return e.getMessage();
        }
    }


    public String getParserDescription() {
        return "Specifies a location with an optional longitude and latitude";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newStringRule(TREE_FILE, false,
                    "A tree file to calculate the average lineage count through time plot."),
            AttributeRule.newIntegerRule(SKIP, true,
                    "The number of trees to skip at start of file."),
            AttributeRule.newIntegerRule(BIN_COUNT, true,
                    "The number of bins to use in lineage count through time plot."),
            AttributeRule.newDoubleRule(MIN_TIME, true,
                    "The min time to compute lineage count from."),
            AttributeRule.newDoubleRule(MAX_TIME, false,
                    "The max time to compute lineage count to."),
    };

    public Class getReturnType() {
        return String.class;
    }
}