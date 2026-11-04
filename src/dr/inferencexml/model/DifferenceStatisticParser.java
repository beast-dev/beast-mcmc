/*
 * DifferenceStatisticParser.java
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

import dr.inference.model.DifferenceStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class DifferenceStatisticParser extends AbstractXMLObjectParser {

    public static String DIFFERENCE_STATISTIC = "differenceStatistic";
    public static String DIFFERENCE = "difference";
    public static String ABSOLUTE = "absolute";

    public String[] getParserNames() {
        return new String[]{getParserName(), DIFFERENCE};
    }

    public String getParserName() {
        return DIFFERENCE_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean absolute = xo.getAttribute(ABSOLUTE, false);

        Statistic term1 = (Statistic) xo.getChild(0);
        Statistic term2 = (Statistic) xo.getChild(1);

        DifferenceStatistic differenceStatistic;
        try {
            differenceStatistic = new DifferenceStatistic(DIFFERENCE_STATISTIC, term1, term2, absolute);
        } catch (IllegalArgumentException iae) {
            throw new XMLParseException("Error parsing " + getParserName() + " the left and right statistics " +
                    "should be of the same dimension or of dimension 1");
        }
        return differenceStatistic;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a statistic that is the difference of the 2 child statistics.";
    }

    public Class getReturnType() {
        return DifferenceStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(ABSOLUTE, true),
            new ElementRule(Statistic.class, "The two operand statistics", 2, 2)
    };

}
