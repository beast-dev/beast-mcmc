/*
 * LoggableStatistic.java
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

package dr.inference.model;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class LoggableStatistic extends Statistic.Abstract {

    public static final String LOGGABLE_STATISTIC = "loggableStatistic";

    private final NumberColumn[] logColumns;

    public LoggableStatistic(String name, NumberColumn[] columns) {
        super(name);
        logColumns = columns;
    }

    public int getDimension() {
        return logColumns.length;
    }

    public double getStatisticValue(int dim) {
        return logColumns[dim].getDoubleValue();
    }

    public String getDimensionName(int dim) {
        return logColumns[dim].getLabel();
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return LOGGABLE_STATISTIC;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String name;
            if (xo.hasAttribute(NAME) || xo.hasAttribute(dr.xml.XMLParser.ID))
                name = xo.getAttribute(NAME, xo.getId());
            else
                name = "";


            final Loggable loggable = (Loggable) xo.getChild(Loggable.class);
            final LogColumn[] logColumns = loggable.getColumns();
            final NumberColumn[] numberColumns = new NumberColumn[logColumns.length];

            for (int i = 0; i < logColumns.length; i++) {
                if (logColumns[i] instanceof NumberColumn) {
                    numberColumns[i] = (NumberColumn) logColumns[i];
                } else
                    throw new XMLParseException("Can only convert NumberColumns into Statistics");
            }

            return new LoggableStatistic(name, numberColumns);
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
                    new ElementRule(Loggable.class),
            };
        }
    };
}
