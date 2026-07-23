/*
 * ReciprocalStatisticParser.java
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

import dr.inference.model.ReciprocalStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class ReciprocalStatisticParser extends AbstractXMLObjectParser {

    public static String RECIPROCAL_STATISTIC = "reciprocalStatistic";
    public static String RECIPROCAL = "reciprocal";

    public String[] getParserNames() {
        return new String[]{getParserName(), RECIPROCAL};
    }

    public String getParserName() {
        return RECIPROCAL_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        ReciprocalStatistic recipStatistic = null;

        Object child = xo.getChild(0);
        if (child instanceof Statistic) {
            recipStatistic = new ReciprocalStatistic((xo.hasId() ? xo.getId() : "RECIPROCAL"), (Statistic) child);
        } else {
            throw new XMLParseException("Unknown element found in " + getParserName() + " element:" + child);
        }

        return recipStatistic;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a statistic that is the element-wise reciprocal of the child statistic.";
    }

    public Class getReturnType() {
        return ReciprocalStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(Statistic.class, 1, 1)
    };
}
