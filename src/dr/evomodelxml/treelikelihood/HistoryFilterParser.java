/*
 * HistoryFilterParser.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.treelikelihood;

import dr.evomodel.treelikelihood.utilities.HistoryFilter;
import dr.util.Identifiable;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class HistoryFilterParser extends AbstractXMLObjectParser {

    public static final String NAME = "historyFilter";
    public static final String MAX_TIME = "maxTime";
    public static final String MIN_TIME = "minTime";
    public static final String SOURCES = "sources";
    public static final String DESTINATIONS = "destinations";
    public static final String INCLUDE_ALL = "includeAll";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double maxTime = xo.getAttribute(MAX_TIME, Double.POSITIVE_INFINITY);
        double minTime = xo.getAttribute(MIN_TIME, 0.0);

        return new HistoryFilter.SetFilter(null, null, maxTime, minTime);
    }

    public String getParserName() {
        return NAME;
    }

    public String getParserDescription() {
        return "A logger to filter transitions in the complete history.";
    }

    public Class getReturnType() {
        return HistoryFilter.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MAX_TIME, true),
            AttributeRule.newDoubleRule(MIN_TIME, true),
            new ElementRule(SOURCES,
                    new XMLSyntaxRule[]{
                            new XORRule(
                                    AttributeRule.newBooleanRule(INCLUDE_ALL, true),
                                    new ElementRule(Identifiable.class, 1, Integer.MAX_VALUE) // TODO Fix type
                            ),
                    }
                    , true),
            new ElementRule(DESTINATIONS,
                    new XMLSyntaxRule[]{
                            new XORRule(
                                    AttributeRule.newBooleanRule(INCLUDE_ALL, true),
                                    new ElementRule(Identifiable.class, 1, Integer.MAX_VALUE) // TODO Fix type
                            ),
                    }
                    , true),
    };
}
