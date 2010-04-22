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

import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

/**
 * A class that stores the distribution statistics for a trace
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TraceDistribution.java,v 1.1.1.2 2006/04/25 23:00:09 rambaut Exp $
 */
public class TraceDistribution<T> {

    public TraceDistribution(T[] values) {
        analyseDistribution(values);
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

        if (values[0] instanceof Double) {
            double[] doubleValues = new double[values.length];
            for (int i = 0; i < values.length; i++) {
                doubleValues[i] = ((Double) values[i]).doubleValue();
            }

            return DiscreteStatistics.meanSquaredError(doubleValues, trueValue);

        } else {
            throw new RuntimeException("Require Continuous Trace Type in the Trace Distribution: " + this);
        }
    }

    private void analyseDistribution(T[] values) {
        this.values = values;

        if (values[0] instanceof Double) {
            double[] doubleValues = new double[values.length];
            for (int i = 0; i < values.length; i++) {
                doubleValues[i] = ((Double) values[i]).doubleValue();
            }
            analyseDistributionContinuous(doubleValues, 0.95);

        } else if (values[0] instanceof Integer) {


        } else if (values[0] instanceof String) {


        } else {
            throw new RuntimeException("Trace type is not recognized: " + values[0].getClass());
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


    public class CredibleSet<T> {
        // <T, frequency> for T = Integer and String
        public Map<T, Integer> credibleSet = new HashMap<T, Integer>();
        public Map<T, Integer> inCredibleSet = new HashMap<T, Integer>();

        public CredibleSet(T[] values, double proportion) {
            credibleSet.clear();
            for (T value : values) {
               if (credibleSet.containsKey(value)) {
                   Integer i = credibleSet.get(value);
                   credibleSet.put(value, i++);
               } else {
                   credibleSet.put(value, 1);
               }
            }

            for (T value : credibleSet.keySet()) {
                double prob = credibleSet.get(value) / values.length;
                if (prob < (1 - proportion)) {
                    inCredibleSet.put(value, credibleSet.get(value));
                    credibleSet.remove(value);
                }
            }
        }

        public CredibleSet(Double[] values, double proportion) {
            analyseDistributionContinuous(Trace.arrayConvert(values), proportion);
        }

        public boolean inside(T value) {
            return credibleSet.containsKey(value);
        }

        public boolean inside(Double value) {
            return value <= hpdUpper && value >= hpdLower;
        }

        private String getSet(Map set) {
            String line = "{";
            Set<T> values = set.keySet();
            TreeSet<T> sortedValues = new TreeSet<T>(values);
            for (T value : sortedValues) {
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