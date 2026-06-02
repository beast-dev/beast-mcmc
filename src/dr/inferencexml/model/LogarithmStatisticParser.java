/*
 * LogarithmStatisticParser.java
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

import dr.inference.model.LogarithmStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class LogarithmStatisticParser extends AbstractXMLObjectParser {

    public static String LOGARITHM_STATISTIC = "logarithmStatistic";
    public static String LOGARITHM = "logarithm";
    public static String BASE = "base";

    public String[] getParserNames() {
        return new String[]{getParserName(), LOGARITHM};
    }

    public String getParserName() {
        return LOGARITHM_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        LogarithmStatistic logStatistic;

        // base 0.0 means natural logarithm, the default base
        double base = xo.getAttribute(BASE, 0.0);

        if (base <= 1.0 && base != 0.0) {
            throw new XMLParseException("Error parsing " + getParserName() + " element: base attribute should be > 1");
        }

        Object child = xo.getChild(0);
        if (child instanceof Statistic) {
            logStatistic = new LogarithmStatistic(LOGARITHM_STATISTIC, (Statistic) child, base);
        } else {
            throw new XMLParseException("Unknown element found in " + getParserName() + " element:" + child);
        }

        return logStatistic;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a statistic that is the element-wise natural logarithm of the child statistic.";
    }

    public Class getReturnType() {
        return LogarithmStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(BASE, true, "An optional base for the logarithm (default is the natural logarithm, base e)"),
            new ElementRule(Statistic.class, 1, 1)
    };

}
