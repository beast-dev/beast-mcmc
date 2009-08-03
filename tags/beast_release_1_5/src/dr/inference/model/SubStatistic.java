/*
 * SubStatistic.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inference.model;

import dr.xml.*;

public class SubStatistic extends Statistic.Abstract {

    public static final String SUB_STATISTIC = "subStatistic";
    public static final String DIMENSION = "dimension";

    private final int[] dimensions;
    private final Statistic statistic;

    public SubStatistic(String name, int[] dimensions, Statistic stat) {
        super(name);
        this.dimensions = dimensions;
        this.statistic = stat;
    }

    public int getDimension() {
        return dimensions.length;
    }

    public double getStatisticValue(int dim) {
        return statistic.getStatisticValue(dimensions[dim]);
    }

    public String getDimensionName(int dim) {
        return statistic.getDimensionName(dimensions[dim]);
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return SUB_STATISTIC;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String name;
            if (xo.hasAttribute(NAME) || xo.hasAttribute(dr.xml.XMLParser.ID))
                name = xo.getAttribute(NAME, xo.getId());
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
    };
}
