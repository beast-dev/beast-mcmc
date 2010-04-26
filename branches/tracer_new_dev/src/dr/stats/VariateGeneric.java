/*
 * Variate.java
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

package dr.stats;

import dr.inference.trace.Trace;

/**
 * interface for a continuous variate.
 *
 * @author Andrew Rambaut
 * @author Walter Xie
 * @version $Id: Variate.java,v 1.15 2006/02/20 17:36:23 rambaut Exp $
 */
public interface VariateGeneric<T> {
    /**
     * @return the number of values in this Variate
     */
    int getCount();

    /**
     * @return the ith value. Returns NaN if it is missing
     */
    T get(int index);

    /**
     * add a value
     */
    void add(T value);

    /**
     * add all the values
     */
    void add(VariateGeneric<T> values);

    /**
     * add all the values
     */
    void add(T[] values);

    /**
     * remove a value
     */
    void remove(int index);

    /**
     * remove all values
     */
    void removeAll();

    /**
     * A String precision concrete implementation of Variate
     */
    public class String implements VariateGeneric<String> {
        private String[] values = null;
        private int valueCount = 0;

        public String() {
        }

        public String(String[] values) {
            this.values = values;
            valueCount = values.length;
        }

        /**
         * add a value
         */
        public void add(String value) {

            if (values == null || valueCount == values.length) {
                String[] newValues = new String[valueCount + 100];
                if (values != null)
                    System.arraycopy(values, 0, newValues, 0, values.length);
                values = newValues;
            }

            values[valueCount] = value;
            valueCount++;
        }

        /**
         * add all the values
         */
        public void add(VariateGeneric<String> values) {
            for (int i = 0; i < values.getCount(); i++) {
                add(values.get(i));
            }
        }

        /**
         * add all the values
         */
        public void add(String[] values) {
            for (String value : values) {
                add(value);
            }
        }


        /**
         * remove a value
         */
        public void remove(int index) {
            if (valueCount > 0) {
                for (int i = index; i < values.length - 1; i++) {
                    values[i] = values[i + 1];
                }
                valueCount--;
            }
        }

        /**
         * remove all values
         */
        public void removeAll() {
            valueCount = 0;
        }

        /**
         * @return the number of values in this Variate
         */
        public int getCount() {
            return valueCount;
        }

        /**
         * @return the ith value.
         */
        public String get(int index) {
            return values[index];
        }


    }

    /**
     * A Integer precision concrete implementation of Variate
     */
    public class Discrete implements VariateGeneric<Integer> {
        private Integer[] values = null;
        private int valueCount = 0;

        public Discrete() {
        }

        public Discrete(Integer[] values) {
            this.values = values;
            valueCount = values.length;
        }

        /**
         * add a value
         */
        public void add(Integer value) {

            if (values == null || valueCount == values.length) {
                Integer[] newValues = new Integer[valueCount + 100];
                if (values != null)
                    System.arraycopy(values, 0, newValues, 0, values.length);
                values = newValues;
            }

            values[valueCount] = value;
            valueCount++;
        }

        /**
         * add all the values
         */
        public void add(VariateGeneric<Integer> values) {
            for (int i = 0; i < values.getCount(); i++) {
                add(values.get(i));
            }
        }

        /**
         * add all the values
         */
        public void add(Integer[] values) {
            for (Integer value : values) {
                add(value);
            }
        }


        /**
         * remove a value
         */
        public void remove(int index) {
            if (valueCount > 0) {
                for (int i = index; i < values.length - 1; i++) {
                    values[i] = values[i + 1];
                }
                valueCount--;
            }
        }

        /**
         * remove all values
         */
        public void removeAll() {
            valueCount = 0;
        }

        /**
         * @return the number of values in this Variate
         */
        public int getCount() {
            return valueCount;
        }

        /**
         * @return the ith value.
         */
        public Integer get(int index) {
            return values[index];
        }

                /**
         * @return the minimum value.
         */
        public Integer getMin() {
            Integer minValue = Integer.MIN_VALUE;

            for (int i = 0; i < valueCount; i++) {
                if (values[i] < minValue)
                    minValue = values[i];
            }

            return minValue;
        }

        /**
         * @return the maximum value.
         */
        public Integer getMax() {
            Integer maxValue = Integer.MAX_VALUE;

            for (int i = 0; i < valueCount; i++) {
                if (values[i] > maxValue)
                    maxValue = values[i];
            }

            return maxValue;
        }
    }

           /**
     * A double precision concrete implementation of Variate
     */
    public class Continuous implements VariateGeneric<Double> {
        private Double[] values = null;
        private int valueCount = 0;

        public Continuous() {
        }

        public Continuous(Double[] values) {
            this.values = values;
            valueCount = values.length;
        }

        /**
         * add a value
         */
        public void add(Double value) {

            if (values == null || valueCount == values.length) {
                Double[] newValues = new Double[valueCount + 100];
                if (values != null)
                    System.arraycopy(values, 0, newValues, 0, values.length);
                values = newValues;
            }

            values[valueCount] = value;
            valueCount++;
        }

        /**
         * add all the values
         */
        public void add(VariateGeneric<Double> values) {
            for (int i = 0; i < values.getCount(); i++) {
                add(values.get(i));
            }
        }

        /**
         * add all the values
         */
        public void add(Double[] values) {
            for (Double value : values) {
                add(value);
            }
        }


        /**
         * remove a value
         */
        public void remove(int index) {
            if (valueCount > 0) {
                for (int i = index; i < values.length - 1; i++) {
                    values[i] = values[i + 1];
                }
                valueCount--;
            }
        }

        /**
         * remove all values
         */
        public void removeAll() {
            valueCount = 0;
        }

        /**
         * @return the number of values in this Variate
         */
        public int getCount() {
            return valueCount;
        }

        /**
         * @return the ith value.
         */
        public Double get(int index) {
            return values[index];
        }

        /**
         * @return if the value for index is missing
         */
        public boolean isMissing(int index) {
            return (values[index].equals(java.lang.Double.NaN));
        }

        /**
         * @return the minimum value.
         */
        public Double getMin() {
            Double minValue = Double.POSITIVE_INFINITY;

            for (int i = 0; i < valueCount; i++) {
                if (values[i] < minValue)
                    minValue = values[i];
            }

            return minValue;
        }

        /**
         * @return the maximum value.
         */
        public Double getMax() {
            Double maxValue = Double.NEGATIVE_INFINITY;

            for (int i = 0; i < valueCount; i++) {
                if (values[i] > maxValue)
                    maxValue = values[i];
            }

            return maxValue;
        }

        /**
         * @return the range.
         */
        public Double getRange() {
            return getMin() - getMax();
        }

        /**
         * @return the mean.
         */
        public Double getMean() {
            return getSum() / getCount();
        }

        /**
         * @return the sum.
         */
        public Double getSum() {
            Double sum = 0.0;

            for (int i = 0; i < getCount(); i++) {
                sum += get(i);
            }

            return sum;
        }

        /**
         * compute the q-th quantile for the distribution
         * (= inverse cdf)
         *
         * @param q quantile (0 <= q <= 1)
         */
        public double getQuantile(Double q) {
            return DiscreteStatistics.quantile(q.doubleValue(), Trace.arrayConvert(values), valueCount);
        }
           }
}
