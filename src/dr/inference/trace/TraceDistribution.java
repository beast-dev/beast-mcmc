/*
 * TraceDistribution.java
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

package dr.inference.trace;

import dr.stats.DiscreteStatistics;
import dr.stats.FrequencyCounter;
import dr.util.HeapSort;

import java.util.*;

/**
 * A class that stores the distribution statistics for a trace
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TraceDistribution.java,v 1.1.1.2 2006/04/25 23:00:09 rambaut Exp $
 */
public class TraceDistribution {
    private TraceType traceType;

    // For discrete and categorical traces
    private Map<Integer, String> categoryLabelMap;
    private List<Integer> categoryOrder;
    private FrequencyCounter<Integer> frequencyCounter;

    public TraceDistribution(List<Double> values, TraceType traceType) {
        this.traceType = traceType;

        assert traceType != TraceType.CATEGORICAL : "cant use this constructor with categorical data";

        initStatistics(values, 0.95);
    }

    public TraceDistribution(List<Double> values, Map<Integer, String> categoryLabelMap, List<Integer> categoryOrder) {
        this.traceType = TraceType.CATEGORICAL;
        this.categoryLabelMap = categoryLabelMap;
        this.categoryOrder = categoryOrder;
        initStatistics(values, 0.95);
    }

    private void initStatistics(List<Double> values, double proportion) {
        if (values.size() < 1) throw new RuntimeException("There is no value sent to statistics calculation !");

        if (traceType.isNumber()) {
            analyseDistributionNumeric(values, proportion);
        }

        if (traceType != TraceType.REAL) {
            analyseDistributionDiscrete(values, proportion);
        }
    }

    public TraceType getTraceType() {
        return traceType;
    }

    public boolean isMinEqualToMax() {
        return minEqualToMax;
    }

    public int getSize() {
        return size;
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

    public double getQ1() {
        return q1;
    }

    public double getQ3() {
        return q3;
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

    public double getMinimum() {
        return minimum;
    }

    public double getMaximum() {
        return maximum;
    }

    public double getHpdLowerCustom() {
        return hpdLowerCustom;
    }

    public double getHpdUpperCustom() {
        return hpdUpperCustom;
    }

    public double getMeanSquaredError(double[] values, double trueValue) {

        if (values == null) {
            throw new RuntimeException("Trace values not yet set");
        }

        if (traceType.isNumber()) {
            return DiscreteStatistics.meanSquaredError(values, trueValue);
        } else {
            throw new RuntimeException("Require Real or Ordinal Trace Type in the Trace Distribution: " + this);
        }
    }

    /**
     * The major method to analyse traces in numeric values including Double, Integer
     * @param valueList the values to analyze
     */
    private void analyseDistributionNumeric(List<Double> valueList, double proportion) {
        double[] values = new double[valueList.size()];
        for (int i = 0; i < valueList.size(); i++) {
            values[i] = valueList.get(i);
        }

        size = values.length;
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
            minEqualToMax = false;
            return;
        }

        int[] indices = new int[values.length];
        HeapSort.sort(values, indices);
        median = DiscreteStatistics.quantile(0.5, values, indices);
        cpdLower = DiscreteStatistics.quantile(0.025, values, indices);
        cpdUpper = DiscreteStatistics.quantile(0.975, values, indices);
        q1 = DiscreteStatistics.quantile(0.25, values, indices);
        q3 = DiscreteStatistics.quantile(0.75, values, indices);
        calculateHPDInterval(proportion, values, indices);
//        ESS = values.length; // move to TraceCorrelation
        calculateHPDIntervalCustom(0.5, values, indices);

        minEqualToMax = true;
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

    private void calculateHPDIntervalCustom(double proportion, double[] array, int[] indices) {
        final double[] hpd = DiscreteStatistics.HPDInterval(proportion, array, indices);
        hpdLowerCustom = hpd[0];
        hpdUpperCustom = hpd[1];
    }

    protected boolean minEqualToMax = false;
    protected boolean hasGeometricMean = false;

    protected int size = 0;
    protected double minimum, maximum;
    protected double mean;
    protected double median, q1, q3;
    protected double geometricMean;
    protected double stdError, meanSquaredError;
    protected double variance;
    protected double cpdLower, cpdUpper, hpdLower, hpdUpper;
    protected double hpdLowerCustom, hpdUpperCustom;

    //************************************************************************
    // new types
    //************************************************************************


    // init FrequencyCounter used for Integer and String
    private void analyseDistributionDiscrete(List<Double> valueList, double proportion) {
        List<Integer> integerValues = new ArrayList<Integer>();
        for (Double value : valueList) {
            integerValues.add(value.intValue());
        }

        if (size == 0) {
            size = integerValues.size();
        }

        if (getTraceType() == TraceType.CATEGORICAL) {
            frequencyCounter = new FrequencyCounter<Integer>(integerValues, proportion);
        } else {
            frequencyCounter = new FrequencyCounter<Integer>(integerValues);
        }

    }

    public Set<Integer> getValueSet() {
        if (categoryOrder != null) {
            return new LinkedHashSet<Integer>(categoryOrder);
        } else {
            return new LinkedHashSet<Integer>(frequencyCounter.getUniqueValues());
        }
    }

    public FrequencyCounter<Integer> getFrequencyCounter() {
        assert traceType.isDiscrete();
        return frequencyCounter;
    }

    public Set<Integer> getCredibleSet() {
        return frequencyCounter.getCredibleSet();
    }

    public Set<Integer> getIncredibleSet() {
        return frequencyCounter.getIncredibleSet();
    }

    public String valueToString(int value) {
        if (categoryLabelMap != null) {
            return categoryLabelMap.get(value);
        }
        return Integer.toString(value);
    }
    public String setToString(Set<Integer> aSet) {
        StringBuilder sb = new StringBuilder("{");
        boolean isFirst = true;
        for (int value : aSet) {
            String label = valueToString(value);
            if (!isFirst) {
                sb.append(", ");
            } else {
                isFirst = false;
            }
            sb.append(label);
        }
        sb.append("}");
        return sb.toString();
    }

    public int getMode() {
        return frequencyCounter.getMode();
    }

    public int getFrequencyOfMode() {
        return frequencyCounter.getFrequency(getMode());
    }

    public double getProbabilityOfMode() {
        if (getSize() > 0) {
            return (double) getFrequencyOfMode() / (double) getSize();
        }
        return 0;
    }

    public List<String> getRange() {
        List<String> valuesList = new ArrayList<String>();
        for (Integer value : frequencyCounter.getUniqueValues() ) {
            if (traceType.isInteger()) { // as Integer is stored as Double in Trace
                if (!valuesList.contains(Integer.toString(((Number) value).intValue())))
                    valuesList.add(Integer.toString(((Number) value).intValue()));
            } else {
                if (!valuesList.contains(value.toString()))
                    valuesList.add(value.toString());
            }
        }
        return valuesList;
    }


}