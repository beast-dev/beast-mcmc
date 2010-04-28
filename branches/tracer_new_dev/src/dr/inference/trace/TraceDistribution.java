/*
 * TraceDistribution.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.inference.trace;

import dr.stats.DiscreteStatistics;
import dr.util.HeapSort;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A class that stores the distribution statistics for a trace
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TraceDistribution.java,v 1.1.1.2 2006/04/25 23:00:09 rambaut Exp $
 */
public class TraceDistribution<T> {

    public TraceDistribution(T[] values) {
        this.values = values;
        credSet = new CredibleSet(values, 0.95);
    }

    public TraceDistribution(T[] values, double ESS) {
        this(values);
        this.ESS = ESS;
    }

    public TraceFactory.TraceType getTraceType() {
        if (values[0].getClass() == TraceFactory.TraceType.CONTINUOUS.getType()) {
            return TraceFactory.TraceType.CONTINUOUS;

        } else if (values[0].getClass() == TraceFactory.TraceType.DISCRETE.getType()) {
            return TraceFactory.TraceType.DISCRETE;

        } else if (values[0].getClass() == TraceFactory.TraceType.CATEGORY.getType()) {
            return TraceFactory.TraceType.CATEGORY;

        } else {
            throw new RuntimeException("Trace type is not recognized: " + values[0].getClass());
        }
    }

    public boolean isValid() {
        return isValid;
    }

    public double getMean() {
        return mean;
    }

    public double getVariance() {
        return variance;
    }

    public double getStdError() {
        return stdError;
    }

    public boolean hasGeometricMean() {
        return hasGeometricMean;
    }

    public double getGeometricMean() {
        return geometricMean;
    }


    public double getMedian() {
        return median;
    }

    public double getLowerHPD() {
        return hpdLower;
    }

    public double getUpperHPD() {
        return hpdUpper;
    }

    public double getLowerCPD() {
        return cpdLower;
    }

    public double getUpperCPD() {
        return cpdUpper;
    }

    public double getESS() {
        return ESS;
    }

    public double getMinimum() {
        return minimum;
    }

    public double getMaximum() {
        return maximum;
    }

    public double getMeanSquaredError(double trueValue) {

        if (values == null) {
            throw new RuntimeException("Trace values not yet set");
        }

        if (values[0] instanceof Number) {
            double[] doubleValues = new double[values.length];
            for (int i = 0; i < values.length; i++) {
                doubleValues[i] = ((Number) values[i]).doubleValue();
            }

            return DiscreteStatistics.meanSquaredError(doubleValues, trueValue);

        } else {
            throw new RuntimeException("Require Number Trace Type in the Trace Distribution: " + this);
        }
    }

    /**
     * @param values the values to analyze
     */
    private void analyseDistributionContinuous(double[] values, double proportion) {
//        this.values = values;   // move to analyseDistribution(T[] values)

        mean = DiscreteStatistics.mean(values);
        stdError = DiscreteStatistics.stdev(values);
        variance = DiscreteStatistics.variance(values);

        minimum = Double.POSITIVE_INFINITY;
        maximum = Double.NEGATIVE_INFINITY;

        for (double value : values) {
            if (value < minimum) minimum = value;
            if (value > maximum) maximum = value;
        }

        if (minimum > 0) {
            geometricMean = DiscreteStatistics.geometricMean(values);
            hasGeometricMean = true;
        }

        if (maximum == minimum) {
            isValid = false;
            return;
        }

        int[] indices = new int[values.length];
        HeapSort.sort(values, indices);
        median = DiscreteStatistics.quantile(0.5, values, indices);
        cpdLower = DiscreteStatistics.quantile(0.025, values, indices);
        cpdUpper = DiscreteStatistics.quantile(0.975, values, indices);
        calculateHPDInterval(proportion, values, indices);
        ESS = values.length;

        isValid = true;
    }

    /**
     * @param proportion the proportion of probability mass included within interval.
     * @param array      the data array
     * @param indices    the indices of the ranks of the values (sort order)
     */
    private void calculateHPDInterval(double proportion, double[] array, int[] indices) {
        final double[] hpd = DiscreteStatistics.HPDInterval(proportion, array, indices);
        hpdLower = hpd[0];
        hpdUpper = hpd[1];
    }

    //************************************************************************
    // private methods
    //************************************************************************

    protected boolean isValid = false;
    protected boolean hasGeometricMean = false;

    protected double minimum, maximum;
    protected double mean;
    protected double median;
    protected double geometricMean;
    protected double stdError;
    protected double variance;
    protected double cpdLower, cpdUpper, hpdLower, hpdUpper;
    protected double ESS;

    protected T[] values;
    public CredibleSet credSet = null;

    public class CredibleSet<T> {
        // <T, frequency> for T = Integer and String
        public Map<T, Integer> valuesMap = new HashMap<T, Integer>();
//        public Map<T, Integer> inCredibleSet = new HashMap<T, Integer>();
        public List<T> credibleSet = new ArrayList<T>();
        public List<T> inCredibleSet = new ArrayList<T>();

        public T mode;
        public int freqOfMode = 0;

        public CredibleSet(T[] values, double proportion) {
            valuesMap.clear();
            credibleSet.clear();
            inCredibleSet.clear();

            if (values[0] instanceof Double) {
                double[] newValues = new double[values.length];
                for (int i = 0; i < values.length; i++) {
                    newValues[i] = ((Double) values[i]).doubleValue();
                }
                analyseDistributionContinuous(newValues, proportion);
            } else {
                if (values[0] instanceof Integer) {
                    double[] newValues = new double[values.length];
                    for (int i = 0; i < values.length; i++) {
                        newValues[i] = ((Integer) values[i]).doubleValue();
                    }
                    analyseDistributionContinuous(newValues, proportion);
                }

                for (T value : values) {
                    if (valuesMap.containsKey(value)) {
                        int i = valuesMap.get(value) + 1;
                        valuesMap.put(value, i);
                    } else {
                        valuesMap.put(value, 1);
                    }
                }

                for (T value : valuesMap.keySet()) {
                    double prob = (double) valuesMap.get(value) / (double) values.length;
                    if (prob < (1 - proportion)) {
                        inCredibleSet.add(value);
                    } else {
                        credibleSet.add(value);
                    }
                }

                calculateMode();
            }
        }

        public boolean inside(T value) {
            return valuesMap.containsKey(value);
        }

        public boolean inside(Double value) {
            return value <= hpdUpper && value >= hpdLower;
        }

        public int getIndex(T value) {
            int i = -1;
            for (T v : valuesMap.keySet()) {
                i++;
                if (v.equals(value)) return i;
            }
            return i;
        }

        public T getMode() {
            return mode;
        }

        public int getFrequencyOfMode() {
            return freqOfMode;
        }

        private void calculateMode() {
            for (T value : valuesMap.keySet()) {
                if (freqOfMode < valuesMap.get(value)) {
                    freqOfMode = valuesMap.get(value);
                    mode = value;
                }
            }
        }

        private String getSet(List<T> list) {
            String line = "{";
            for (T value : list) {
                line = line + value + ", ";
            }
            if (line.endsWith(", ")) {
                line = line.substring(0, line.lastIndexOf(", ")) + "}";
            } else {
                line = "n/a";
            }
            return line;
        }

        public String getCredibleSet() {
            return getSet(credibleSet);
        }

        public String getInCredibleSet() {
            return getSet(inCredibleSet);
        }
    }


}