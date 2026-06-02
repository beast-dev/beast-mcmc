/*
 * SubStatisticParser.java
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

import dr.inference.model.Statistic;
import dr.inference.model.SubStatistic;
import dr.xml.*;

/**
 */
public class SubStatisticParser extends AbstractXMLObjectParser {

    public static final String SUB_STATISTIC = "subStatistic";
    public static final String DIMENSION = "dimension";

    public String getParserName() {
        return SUB_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name;
        if (xo.hasAttribute(Statistic.NAME) || xo.hasAttribute(dr.xml.XMLParser.ID))
            name = xo.getAttribute(Statistic.NAME, xo.getId());
        else
            name = "";

        final Statistic stat = (Statistic) xo.getChild(Statistic.class);

        final int[] values = xo.getIntegerArrayAttribute(DIMENSION);

        if (values.length == 0) {
            throw new XMLParseException("Must specify at least one dimension");
        }

        final int dim = stat.getDimension();

        for (int value : values) {
            if (value >= dim || value < 0) {
                throw new XMLParseException("Dimension " + value + " is not a valid dimension.");
            }
        }

        return new SubStatistic(name, values, stat);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Allows you to choose specific dimensions of a given statistic";
    }

    public Class getReturnType() {
        return SubStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newIntegerArrayRule(DIMENSION, false),
                new ElementRule(Statistic.class),
        };
    }
}
