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

    private Map<Integer, String> categoryLabelMap;
    private Map<Integer, Integer> categoryOrderMap;

    public TraceDistribution(List<Double> values, TraceType traceType) {
        this.traceType = traceType;

        assert traceType != TraceType.CATEGORICAL : "cant use this constructor with categorical data";

        initStatistics(values, 0.95);
    }

    public TraceDistribution(List<Double> values, Map<Integer, String> categoryLabelMap) {
        this.traceType = TraceType.CATEGORICAL;
        this.categoryLabelMap = categoryLabelMap;
        this.categoryOrderMap = getNaturalOrder(categoryLabelMap);
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


    public Map<Integer, Integer> getNaturalOrder(Map<Integer, String> categoryMap) {
        categoryOrderMap = new HashMap<Integer, Integer>();

        List<String> labels = new ArrayList<String>(categoryMap.values());
        Collections.sort(labels);

        for (Integer index : categoryMap.keySet()) {
            String l = categoryMap.get(index);
            categoryOrderMap.put(labels.indexOf(l), index);
        }

        return categoryOrderMap;
    }

    public Map<Integer, Integer> getFrequencyOrder() {
        categoryOrderMap = new HashMap<Integer, Integer>();

        List<Integer> values = frequencyCounter.getUniqueValues();
        return categoryOrderMap;
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

    // frequency counter for discrete traces
    private FrequencyCounter<Integer> frequencyCounter;

    // init FrequencyCounter used for Integer and String
    private void analyseDistributionDiscrete(List<Double> valueList, double proportion) {
        List<Integer> values = new ArrayList<Integer>();
        for (Double value : valueList) {
            values.add(value.intValue());
        }

        if (size == 0) {
            size = values.size();
        }
        frequencyCounter = new FrequencyCounter<Integer>(values);
    }

    // categorical data is stored as integer values (in a double list) - mapping to
    // category labels is done by the Traces themselves.
//    /**
//     * Convert a categorical value to the index of unique values,
//     * and put it into a map <code>Map<Integer, String></code>.
//     * Key is the index of unique values, value of map is that unique value.
//     *
//     * @return
//     */
//    public Map<Integer, String> getIndexMap() {
//        Map<Integer, String> categoryDataMap = new HashMap<Integer, String>();
//        if (frequencies != null && categoryDataMap.size() == 0) {
//            int i = -1;
//            for (Object key : frequencies.uniqueValues()) {
//                i++;
//                String value = key.toString();
//                categoryDataMap.put(i, value);
//            }
//        }
//        return categoryDataMap;
//    }

//    /**
//     * Convert a list of categorical values into a list of indices of their unique values,
//     * given the map <code>categoryDataMap</code>.
//     * Key is the index of unique values, value of map is that unique value
//     *
//     * @param values
//     * @return
//     */
//    public List<Double> indexingData(List<String> values) {
//        List<Double> intData = new ArrayList<Double>();
//        Map<Integer, String> categoryDataMap = getIndexMap();
//        if (categoryDataMap.size() < 1) return intData;
//
//        for (int v = 0; v < values.size(); v++) {
//            for (Map.Entry<Integer, String> entry : categoryDataMap.entrySet()) {
//                if (values.get(v).equals(entry.getValue())) {
//                    // add index
//                    intData.add(v, (double) entry.getKey());
//                    break;
//                }
//            }
//
//        }
//        if (intData.size() > 0 && values.size() != intData.size())
//            System.err.println("values.size(" + values.size() + ") != intData.size(" + intData.size() + ") !");
//
//        return intData;
//    }

    public Set<Integer> getValueSet() {
        Set<Integer> valueSet = new LinkedHashSet<Integer>();
        for (Integer value : frequencyCounter.getUniqueValues()) {
            int index = (categoryOrderMap != null ? categoryOrderMap.get(value) : value);
            valueSet.add(index);
        }
        return valueSet;
    }

    public Set<Integer> getCredibleSet() {
        Set<Integer> credibleSet = new LinkedHashSet<Integer>();
        for (Integer value : frequencyCounter.getCredibleSet()) {
            int index = (categoryOrderMap != null ? categoryOrderMap.get(value) : value);
            credibleSet.add(index);
        }
        return credibleSet;
    }

    public Set<Integer> getIncredibleSet() {
        Set<Integer> incredibleSet = new LinkedHashSet<Integer>();
        for (Integer value : frequencyCounter.getIncredibleSet()) {
            int index = (categoryOrderMap != null ? categoryOrderMap.get(value) : value);
            incredibleSet.add(index);
        }
        return incredibleSet;
    }

    public String setToString(Set<Integer> aSet) {
        StringBuilder sb = new StringBuilder("{");
        boolean isFirst = true;
        for (int value : aSet) {
            String label = categoryLabelMap.get(value);
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

//    private boolean contains(Set<T> aSet, int valueORIndex) {
//        if (traceType.isNumber()) {
//            // T is either Double or Integer
//            return aSet.contains((double) valueORIndex);
//        } else { // String
//            String valueString = null;
//            int i = -1;
//            for (T v : frequencies.uniqueValues()) {
//                i++;
//                if (i == valueORIndex) {
//                    valueString = v.toString();
//                    break;
//                }
//            }
//            if (valueString == null) {
//                return false;
//            }
//            return aSet.contains(valueString);
//        }
//    }

}