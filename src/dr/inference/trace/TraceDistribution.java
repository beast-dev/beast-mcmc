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

import dr.stats.CredibleSet;
import dr.stats.DiscreteStatistics;
import dr.stats.FrequencyCounter;
import dr.stats.Mode;
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
    private TraceType traceType;

    public TraceDistribution(List<T> values, TraceType traceType) {
        this.traceType = traceType;
        initStatistics(values, 0.95);
    }

    public TraceDistribution(List<T> values, TraceType traceType, double ESS) {
        this(values, traceType);
        this.ESS = ESS;
    }

    public TraceType getTraceType() {
        return traceType;
    }

    public void setTraceType(TraceType traceType) {
        this.traceType = traceType;
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
     * @param values the values to analyze
     */
    private void analyseDistributionNumeric(double[] values, double proportion) {
//        this.values = values;   // move to TraceDistribution(T[] values)

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
        calculateHPDIntervalCustom(0.5, values, indices);

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

    private void calculateHPDIntervalCustom(double proportion, double[] array, int[] indices) {
        final double[] hpd = DiscreteStatistics.HPDInterval(proportion, array, indices);
        hpdLowerCustom = hpd[0];
        hpdUpperCustom = hpd[1];
    }

    protected boolean isValid = false;
    protected boolean hasGeometricMean = false;

    protected double minimum, maximum;
    protected double mean;
    protected double median;
    protected double geometricMean;
    protected double stdError, meanSquaredError;
    protected double variance;
    protected double cpdLower, cpdUpper, hpdLower, hpdUpper;
    protected double hpdLowerCustom, hpdUpperCustom;
    protected double ESS;

    //************************************************************************
    // new types
    //************************************************************************

    // frequency counter for T = Integer and String
    public FrequencyCounter<T> frequencyCounter;
    protected Mode<T> mode;
    public CredibleSet<T> credibleSet;

    public void initStatistics(List<T> values, double proportion) {
        if (values.size() < 1) throw new RuntimeException("There is no value sent to statistics calculation !");

        if (traceType.isNumber()) {
            double[] newValues = new double[values.size()];
            for (int i = 0; i < values.size(); i++) {
                newValues[i] = ((Number) values.get(i)).doubleValue();
            }
            analyseDistributionNumeric(newValues, proportion);
        }

        if (traceType != TraceType.REAL) {
            analyseDistributionDiscrete(values, proportion);
        }
    }

    // init FrequencyCounter used for Integer and String
    private void analyseDistributionDiscrete(List<T> values, double proportion) {
        frequencyCounter = new FrequencyCounter<T>(values);
        mode = frequencyCounter.getMode();
        credibleSet = frequencyCounter.getCredibleSet(proportion);

        isValid = true; // todo what purpose?
    }

    //************ Used by panels or FrequencyPlot *************
    public int getIndex(T value) {
        return frequencyCounter.getKeyIndex(value);
    }

    public boolean credibleSetContains(int valueORIndex) {
        return contains(credibleSet.getCredibleSet(), valueORIndex);
    }

    public boolean incredibleSetContains(int valueORIndex) {
        return contains(credibleSet.getIncredibleSet(), valueORIndex);
    }

    public String printCredibleSet() {
        return credibleSet.toStringCredibleSet();
    }

    public String printIncredibleSet() {
        return credibleSet.toStringIncredibleSet();
    }

    public T getMode() {
        if (mode == null) return null;
        return mode.getMode();
    }

    public int getFrequencyOfMode() {
        if (mode == null) return 0;
        return mode.getFrequencyOfMode();
    }

    public List<String> getRange() {
        List<String> valuesList = new ArrayList<String>();
        for (T value : frequencyCounter.uniqueValues()) {
            if (traceType == TraceType.ORDINAL) { // as Integer is stored as Double in Trace
                if (!valuesList.contains(Integer.toString(((Number) value).intValue())))
                    valuesList.add(Integer.toString(((Number) value).intValue()));
            } else {
                if (!valuesList.contains(value.toString()))
                    valuesList.add(value.toString());
            }
        }
        return valuesList;
    }

    private boolean contains(Set<T> aSet, int valueORIndex) {
        if (traceType.isNumber()) {
            // T is either Double or Integer
            return aSet.contains((double) valueORIndex);
        } else { // String
            String valueString = null;
            int i = -1;
            for (T v : frequencyCounter.uniqueValues()) {
                i++;
                if (i == valueORIndex) valueString = v.toString();
            }
            return aSet.contains(valueString);
        }
    }

}