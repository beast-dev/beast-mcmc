/*
 * SumStatisticParser.java
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

package dr.inferencexml.model;

import dr.inference.model.Statistic;
import dr.inference.model.SumStatistic;
import dr.xml.*;

/**
 */
public class SumStatisticParser extends AbstractXMLObjectParser {

    public static String SUM_STATISTIC = "sumStatistic";
    public static String SUM = "sum";
    public static String ELEMENTWISE = "elementwise";
    public static String ABSOLUTE = "absolute";

    public String[] getParserNames() {
        return new String[]{getParserName(), SUM};
    }

    public String getParserName() {
        return SUM_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean elementwise = xo.getAttribute(ELEMENTWISE, false);
        boolean absolute = xo.getAttribute(ABSOLUTE, false);

        String name = SUM_STATISTIC;
        if (xo.hasAttribute(Statistic.NAME)) {
            name = xo.getAttribute(Statistic.NAME, xo.getId());
        } else if (xo.hasAttribute(XMLParser.ID)) {
            name = xo.getAttribute(XMLParser.ID, xo.getId());
        }

        final SumStatistic sumStatistic = new SumStatistic(name, elementwise, absolute);

        for (int i = 0; i < xo.getChildCount(); i++) {
            final Statistic statistic = (Statistic) xo.getChild(i);

            try {
                sumStatistic.addStatistic(statistic);
            } catch (IllegalArgumentException iae) {
                throw new XMLParseException("Statistic added to " + getParserName() + " element is not of the same dimension");
            }
        }

        return sumStatistic;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a statistic that is the element-wise sum of the child statistics.";
    }

    public Class getReturnType() {
        return SumStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(ELEMENTWISE, true),
            AttributeRule.newStringRule(Statistic.NAME, true),
            new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
    };
}
