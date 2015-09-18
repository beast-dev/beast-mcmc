/*
 * MultiEpochExponential.java
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

package dr.evolution.coalescent;

/**
 * This class models a multi-phase exponential growth
 *
 * @author Marc A. Suchard
 */

public class MultiEpochExponential extends ConstantPopulation {

    /**
     * Construct demographic model with default settings
     */
    public MultiEpochExponential(Type units, int numEpoch) {

        super(units);
        transitionTime = new double[numEpoch - 1];
        rate = new double[numEpoch];
    }

    public void setTransitionTime(int index, double transitionTime) {
        this.transitionTime[index] = transitionTime;
    }

    public void setGrowthRate(int index, double rate) {
        this.rate[index] = rate;
    }

    public double getDemographic(double t) {

        double logDemographic = Math.log(getN0()); // TODO work on log-scale
        double lastTransitionTime = 0.0;
        int currentEpoch = 0;

        // Account for all epoch until t
        while (currentEpoch < transitionTime.length && t > transitionTime[currentEpoch]) {
            logDemographic += -rate[currentEpoch] * (transitionTime[currentEpoch] - lastTransitionTime);
            lastTransitionTime = transitionTime[currentEpoch];
            ++currentEpoch;
        }

        // Account for epoch that holds t
        logDemographic += -rate[currentEpoch] * (t - lastTransitionTime);

        return Math.exp(logDemographic);
    }

    /**
     * Calculates the integral 1/N(x) dx between start and finish.
     */
    public double getIntegral(double start, double finish) {
        // TODO An analytic solution exists
        return getNumericalIntegral(start, finish);
    }

    public double getIntensity(double t) { throw new RuntimeException("Not implemented!"); }

    public double getInverseIntensity(double x) {
        throw new RuntimeException("Not implemented!");
    }

    public int getNumArguments() {
        throw new RuntimeException("Not implemented!");
    }

    public String getArgumentName(int n) {
        throw new RuntimeException("Not implemented!");
    }

    public double getArgument(int n) {
        throw new RuntimeException("Not implemented!");
    }

    public void setArgument(int n, double value) {
        throw new RuntimeException("Not implemented!");
    }

    public double getLowerBound(int n) {
        throw new RuntimeException("Not implemented!");
    }

    public double getUpperBound(int n) {
        throw new RuntimeException("Not implemented!");
    }

    //
    // private stuff
    //

    final private double[] transitionTime;
    final private double[] rate;
}
