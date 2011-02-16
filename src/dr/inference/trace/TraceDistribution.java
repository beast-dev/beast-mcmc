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

import java.util.*;

/**
 * A class that stores the distribution statistics for a trace
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TraceDistribution.java,v 1.1.1.2 2006/04/25 23:00:09 rambaut Exp $
 */
public class TraceDistribution<T> {

    Class traceType;

    public TraceDistribution(List<T> values) {
        credSet = new CredibleSet(values, 0.95);
    }

    public TraceDistribution(List<T> values, double ESS) {
        this(values);
        this.ESS = ESS;
    }

    public Class getTraceType() {
        return traceType;
    }

    public String getTraceTypeBrief() {
        if (traceType == Double.class) {
            return TraceFactory.TraceType.CONTINUOUS.getBrief();
        } else if (traceType == Integer.class) {
            return TraceFactory.TraceType.INTEGER.getBrief();
        } else if (traceType == String.class) {
            return TraceFactory.TraceType.CATEGORY.getBrief();
        }
        throw new IllegalArgumentException("The trace type " + traceType + " is not recognized.");
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

    public double getMeanSquaredError(double[] values, double trueValue) {

        if (values == null) {
            throw new RuntimeException("Trace values not yet set");
        }

        if (traceType == Number.class) {
            return DiscreteStatistics.meanSquaredError(values, trueValue);
        } else {
            throw new RuntimeException("Require Number Trace Type in the Trace Distribution: " + this);
        }
    }

    /**
     * @param valuesC the values to analyze
     */
    private void analyseDistributionContinuous(double[] valuesC, double proportion) {
//        this.values = values;   // move to TraceDistribution(T[] values)

        mean = DiscreteStatistics.mean(valuesC);
        stdError = DiscreteStatistics.stdev(valuesC);
        variance = DiscreteStatistics.variance(valuesC);

        minimum = Double.POSITIVE_INFINITY;
        maximum = Double.NEGATIVE_INFINITY;

        for (double value : valuesC) {
            if (value < minimum) minimum = value;
            if (value > maximum) maximum = value;
        }

        if (minimum > 0) {
            geometricMean = DiscreteStatistics.geometricMean(valuesC);
            hasGeometricMean = true;
        }

        if (maximum == minimum) {
            isValid = false;
            return;
        }

        int[] indices = new int[valuesC.length];
        HeapSort.sort(valuesC, indices);
        median = DiscreteStatistics.quantile(0.5, valuesC, indices);
        cpdLower = DiscreteStatistics.quantile(0.025, valuesC, indices);
        cpdUpper = DiscreteStatistics.quantile(0.975, valuesC, indices);
        calculateHPDInterval(proportion, valuesC, indices);
        ESS = valuesC.length;

//        isValid = true;
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
    protected double stdError, meanSquaredError;
    protected double variance;
    protected double cpdLower, cpdUpper, hpdLower, hpdUpper;
    protected double ESS;

//    protected T[] values;
    public CredibleSet credSet = null;

    public class CredibleSet<T> {
        // <T, frequency> for T = Integer and String
        public Map<T, Integer> valuesMap = new HashMap<T, Integer>();
//        public Map<T, Integer> inCredibleSet = new HashMap<T, Integer>();
        public List<T> credibleSet = new ArrayList<T>();
        public List<T> inCredibleSet = new ArrayList<T>();

        public T mode;
        public int freqOfMode = 0;

        public CredibleSet(List<T> values, double proportion) {
            valuesMap.clear();
            credibleSet.clear();
            inCredibleSet.clear();

            if (!(values.get(0) instanceof Double)) {// make sure: if T is Object then default to double
                if (values.get(0) instanceof Integer) {
                    traceType = Integer.class;

                    double[] newValues = new double[values.size()];
                    for (int i = 0; i < values.size(); i++) {
                        newValues[i] = ((Integer) values.get(i)).doubleValue();
                    }
                    analyseDistributionContinuous(newValues, proportion);
                } else {
                    traceType = String.class;
                }

                for (T value : values) {
                    if (valuesMap.containsKey(value)) {
                        int i = valuesMap.get(value) + 1;
                        valuesMap.put(value, i);
                    } else {
                        valuesMap.put(value, 1);
                    }
                }

                for (T value : new TreeSet<T>(valuesMap.keySet())) {
                    double prob = valuesMap.get(value).doubleValue() / (double) values.size();
                    if (prob < (1 - proportion)) {
                        inCredibleSet.add(value);
                    } else {
                        credibleSet.add(value);
                    }
                }

                calculateMode();

            } else {
                traceType = Double.class;
                double[] newValues = new double[values.size()];
                for (int i = 0; i < values.size(); i++) {
                    newValues[i] = ((Double) values.get(i)).doubleValue();
                }
                analyseDistributionContinuous(newValues, proportion);
            }

            isValid = true;

        }

        public boolean inside(T value) {
            return valuesMap.containsKey(value);
        }

        public boolean inside(Double value) {
            return value <= hpdUpper && value >= hpdLower;
        }

        public int getIndex(T value) {
            int i = -1;
            for (T v : new TreeSet<T>(valuesMap.keySet())) {
                i++;
                if (v.equals(value)) return i;
            }
            return i;
        }

        public boolean credibleSetContains(int valueORIndex) {
            return contains(credibleSet, valueORIndex);
        }

        public boolean inCredibleSetContains(int valueORIndex) {
            return contains(inCredibleSet, valueORIndex);
        }

        private boolean contains(List<T> list, int valueORIndex) {
            if (traceType == Integer.class) {
                return list.contains(valueORIndex);
            } else { // String
                String valueString = null;
                int i = -1;
                for (T v : new TreeSet<T>(valuesMap.keySet())) {
                    i++;
                    if (i == valueORIndex) valueString = v.toString();
                }
                return list.contains(valueString);
            }
        }

        public T getMode() {
            return mode;
        }

        public int getFrequencyOfMode() {
            return freqOfMode;
        }

        public List<String> getRange() {
            List<String> valuesList = new ArrayList<String>();
            for (T value : new TreeSet<T>(valuesMap.keySet())) {
                if (!valuesList.contains(value.toString()))
                    valuesList.add(value.toString());
            }
            return valuesList;
        }

        private void calculateMode() {
            for (T value : new TreeSet<T>(valuesMap.keySet())) {
                if (freqOfMode < valuesMap.get(value)) {
                    freqOfMode = valuesMap.get(value);
                    mode = value;
                }
            }
        }

        private String printSet(List<T> list) {
            String line = "{";
            for (T value : list) {
                line = line + value + ", ";
            }
            if (line.endsWith(", ")) {
                line = line.substring(0, line.lastIndexOf(", ")) + "}";
            } else {
                line = "{}";
            }
            return line;
        }

        public String printCredibleSet() {
            return printSet(credibleSet);
        }

        public String printInCredibleSet() {
            return printSet(inCredibleSet);
        }
    }
}