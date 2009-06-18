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

/**
 * A class that stores the distribution statistics for a trace
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TraceDistribution.java,v 1.1.1.2 2006/04/25 23:00:09 rambaut Exp $
 */
public class TraceDistribution {

    public TraceDistribution(double[] values) {
        analyseDistribution(values);
    }

    public TraceDistribution(double[] values, double ESS) {
        analyseDistribution(values);
        this.ESS = ESS;
    }

    public boolean isValid() {
        return isValid;
    }

    public double getMean() {
        return mean;
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

    /**
     * @param values the values to analyze
     */
    private void analyseDistribution(double[] values) {

        mean = DiscreteStatistics.mean(values);
        stdError = DiscreteStatistics.stdev(values);

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
        calculateHPDInterval(0.95, values, indices);
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
    protected double mean, median, geometricMean, stdError;
    protected double cpdLower, cpdUpper, hpdLower, hpdUpper;
    protected double ESS;
}