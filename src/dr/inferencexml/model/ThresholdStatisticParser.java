/*
 * ThresholdStatisticParser.java
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

package dr.inferencexml.model;

import dr.inference.model.Statistic;
import dr.inference.model.SumStatistic;
import dr.inference.model.ThresholdStatistic;
import dr.xml.*;

/**
 */
public class ThresholdStatisticParser extends AbstractXMLObjectParser {

    public static final String THRESHOLD_STATISTIC = "thresholdStatistic";

    public static String THRESHOLD = "threshold";

    public String getParserName() {
        return THRESHOLD_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double threshold = xo.getAttribute(THRESHOLD, 0.0);

        String name = THRESHOLD_STATISTIC;
        if (xo.hasAttribute(Statistic.NAME)) {
            name = xo.getAttribute(Statistic.NAME, xo.getId());
        } else if (xo.hasAttribute(XMLParser.ID)) {
            name = xo.getAttribute(XMLParser.ID, xo.getId());
        }

        final Statistic statistic = (Statistic) xo.getChild(Statistic.class);

        final ThresholdStatistic thresholdStatistic = new ThresholdStatistic(name, statistic, threshold);

        return thresholdStatistic;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a statistic that thresholds the child statistics.";
    }

    public Class getReturnType() {
        return ThresholdStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(Statistic.NAME, true),
            new ElementRule(Statistic.class),
            AttributeRule.newDoubleRule(THRESHOLD, true),
    };
}
