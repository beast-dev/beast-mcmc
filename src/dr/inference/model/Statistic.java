/*
 * Statistic.java
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
import dr.util.Attribute;
import dr.util.Identifiable;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: Statistic.java,v 1.8 2005/05/24 20:26:00 rambaut Exp $
 */
public interface Statistic extends Attribute<double[]>, Identifiable, Loggable {

    public static final String NAME = "name";

    /**
     * @return the name of this statistic
     */
    String getStatisticName();

    /**
     * @param dim the dimension to return name of
     * @return the statistic's name for a given dimension
     */
    String getDimensionName(int dim);

    /**
     * Set the names of the dimensions (optional, by default they are named after the statistic).
     * @param names
     */
    void setDimensionNames(String[] names) ;

    /**
     * @return the number of dimensions that this statistic has.
     */
    int getDimension();

    /**
     * @param dim the dimension to return value of
     * @return the statistic's scalar value in the given dimension
     */
    double getStatisticValue(int dim);

    /**
     * @return the sum of all the elements
     */
    double getValueSum();

    /**
     * Abstract base class for Statistics
     */
    public abstract class Abstract implements Statistic {

        private String name = null;

        public Abstract() {
            this.name = null;
        }

        public Abstract(String name) {
            this.name = name;
        }

        @Override
        public String getStatisticName() {
            if (name != null) {
                return name;
            } else if (id != null) {
                return id;
            } else {
                return getClass().toString();
            }
        }

        @Override
        public String getDimensionName(int dim) {
            if (getDimension() == 1) {
                return getStatisticName();
            } else {
                return getStatisticName() + Integer.toString(dim + 1);
            }
        }

        @Override
        public void setDimensionNames(String[] names) {
            // do nothing
        }

        @Override
        public double getValueSum() {
            double sum = 0.0;
            for (int i = 0; i < getDimension(); i++) {
                sum += getStatisticValue(i);
            }
            return sum;
        }

        @Override
        public String toString() {
            StringBuffer buffer = new StringBuffer(String.valueOf(getStatisticValue(0)));

            for (int i = 1; i < getDimension(); i++) {
                buffer.append(", ").append(String.valueOf(getStatisticValue(i)));
            }
            return buffer.toString();
        }

        // **************************************************************
        // Attribute IMPLEMENTATION
        // **************************************************************

        @Override
        public final String getAttributeName() {
            return getStatisticName();
        }

        @Override
        public double[] getAttributeValue() {
            double[] stats = new double[getDimension()];
            for (int i = 0; i < stats.length; i++) {
                stats[i] = getStatisticValue(i);
            }

            return stats;
        }

        // **************************************************************
        // Identifiable IMPLEMENTATION
        // **************************************************************

        protected String id = null;

        /**
         * @return the id.
         */
        @Override
        public String getId() {
            return id;
        }

        /**
         * Sets the id.
         */
        @Override
        public void setId(String id) {
            this.id = id;
        }

        // **************************************************************
        // Loggable IMPLEMENTATION
        // **************************************************************

        /**
         * @return the log columns.
         */
        @Override
        public LogColumn[] getColumns() {
            LogColumn[] columns = new LogColumn[getDimension()];
            for (int i = 0; i < getDimension(); i++) {
                columns[i] = new StatisticColumn(getDimensionName(i), i);
            }
            return columns;
        }

        private class StatisticColumn extends NumberColumn {
            private final int dim;

            public StatisticColumn(String label, int dim) {
                super(label);
                this.dim = dim;
            }

            public double getDoubleValue() {
                return getStatisticValue(dim); }
        }
    }
}
