/*
 * MadStatisticParser.java
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

package dr.inferencexml.model;

import dr.inference.model.MadStatistic;
import dr.inference.model.Statistic;
import dr.inference.model.VarianceStatistic;
import dr.xml.*;

/**
 */
public class MadStatisticParser extends AbstractXMLObjectParser {

    public static String MAD_STATISTIC = "madStatistic";
    public static String CIM = "centerIsMedian";
    public static String DIM = "dispersionIsMedian";

    public String getParserName() { return MAD_STATISTIC; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean centerIsMedian = xo.getAttribute(CIM,true);
        boolean dispersionIsMedian = xo.getAttribute(DIM,true);

        MadStatistic madStatistic = new MadStatistic(MAD_STATISTIC,centerIsMedian,dispersionIsMedian);

        for (int i =0; i < xo.getChildCount(); i++) {
            Object child = xo.getChild(i);
            if (child instanceof Statistic) {
                madStatistic.addStatistic((Statistic)child);
            } else {
                throw new XMLParseException("Unknown element found in " + getParserName() + " element:" + child);
            }
        }

        return madStatistic;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a statistic that is the MAD of the child statistics. Default is Median Absolute Deviation, median(abs(x - median(x))).";
    }

    public Class getReturnType() { return VarianceStatistic.class; }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            new ElementRule(Statistic.class, 1, Integer.MAX_VALUE ),
            AttributeRule.newBooleanRule(CIM, true, "If true, the central tendency taken is center(x) = median(x). If false, center(x) = mean(x). Default true for median."),
            AttributeRule.newBooleanRule(DIM, true, "If true, the dispersion tendency taken is median(abs(x - center(x))). If false, mean(abs(x - center(x))). Default true for median.")
    };

}
