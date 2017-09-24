/*
 * Variate.java
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

package dr.stats;

import java.util.ArrayList;
import java.util.List;

/**
 * interface for a continuous variate.
 *
 * @author Andrew Rambaut
 * @version $Id: Variate.java,v 1.15 2006/02/20 17:36:23 rambaut Exp $
 */
public interface Variate<T> {
    /**
     * @return the number of values in this Variate
     */
    int getCount();

    /**
     * @return the ith value. Returns NaN if it is missing
     */
    T get(int index);

    /**
     * @return the minimum value.
     */
    T getMin();

    /**
     * @return the maximum value.
     */
    T getMax();

    /**
     * @return the range.
     */
    T getRange();

    /**
     * @return the sum.
     */
    T getSum();

    /**
     * @return the mean.
     */
    T getMean();

    /**
     * compute the q-th quantile for the distribution
     * (= inverse cdf)
     *
     * @param q quantile (0 <= q <= 1)
     */
    double getQuantile(T q);

    /**
     * add a value
     */
    void add(T value);

    /**
     * add all the values
     */
    void add(Variate values);

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


     abstract class N<Number> implements Variate<Number> {
        protected List<Number> values = new ArrayList<Number>();

//        public N() {}
//        public N(List<Number> values) {
//            add(values);
//        }

        public int getCount() {
            return values.size();
        }

        public void add(Variate values) {
            for (int i = 0; i < values.getCount(); i++) {
                add((Number) values.get(i));
            }
        }

        public void add(Number value) {
            values.add(value);
        }

        public Number get(int index) {
            return values.get(index);
        }

        public void remove(int index) {
            values.remove(index);
        }

        public void removeAll() {
            values.clear();
        }

        public void add(List<Number> values) {
            this.values.addAll(values);
        }

        public void add(Number[] values) {
            for (Number value : values) {
                add(value);
            }
        }
    }

    /**
     * A double precision concrete implementation of Variate
     */
     class D extends N<Double> {
        public D() {}

        public D(List<Double> values) {
//            super(values);
            add(values);
        }

        public D(Double[] values) {
            add(values);
        }

        /**
         * @return the minimum value.
         */
        public Double getMin() {
            Double minValue = java.lang.Double.POSITIVE_INFINITY;
            for (Double value : values) {
                if (value < minValue)
                    minValue = value;
            }
            return minValue;
        }

        /**
         * @return the maximum value.
         */
        public Double getMax() {
            Double maxValue = java.lang.Double.NEGATIVE_INFINITY;
            for (Double value : values) {
                if (value > maxValue)
                    maxValue = value;
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
            for (Double value : values) {
                sum += value;
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
            double[] dv = new double[values.size()];
            for (int i = 0; i < values.size(); i++) {
                dv[i] = values.get(i);
            }
            return DiscreteStatistics.quantile(q.doubleValue(), dv);
        }

    }

    /**
     * A single precision concrete implementation of Variate
     */
     class I extends N<Integer> {
        public I() {}

        public I(List<Integer> values) {
//            super(values);
            add(values);
        }

        public I(Integer[] values) {
            add(values);
        }

        /**
         * @return the minimum value.
         */
        public Integer getMin() {
            Integer minValue = java.lang.Integer.MAX_VALUE;
            for (Number value : values) {
                if (value.intValue() < minValue)
                    minValue = value.intValue();
            }
            return minValue;
        }

        /**
         * @return the maximum value.
         */
        public Integer getMax() {
            Integer maxValue = java.lang.Integer.MIN_VALUE;
            for (Number value : values) {
                if (value.intValue() > maxValue)
                    maxValue = value.intValue();
            }
            return maxValue;
        }

        /**
         * @return the range.
         */
        public Integer getRange() {
            return getMin() - getMax();
        }

        /**
         * @return the mean.
         */
        public Integer getMean() {
            return getSum() / getCount();
        }

        /**
         * @return the sum.
         */
        public Integer getSum() {
            Integer sum = 0;
            for (Integer value : values) {
                sum += value;
            }
            return sum;
        }

        /**
         * compute the q-th quantile for the distribution
         * (= inverse cdf)
         *
         * @param q quantile (0 <= q <= 1)
         */
        public double getQuantile(Integer q) {
            double[] dv = new double[values.size()];
            for (int i = 0; i < values.size(); i++) {
                dv[i] = values.get(i).doubleValue();
            }
            return DiscreteStatistics.quantile(q.doubleValue(), dv);
        }

    }
}