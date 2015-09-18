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

        double logDemographic = 0.0;
        double lastTransitionTime = 0.0;
        int currentEpoch = 0;

        // Account for all epochs before t
        while (currentEpoch < transitionTime.length && t > transitionTime[currentEpoch]) {
            logDemographic += -rate[currentEpoch] * (transitionTime[currentEpoch] - lastTransitionTime);
            lastTransitionTime = transitionTime[currentEpoch];
            ++currentEpoch;
        }

        // Account for epoch that holds t
        logDemographic += -rate[currentEpoch] * (t - lastTransitionTime);

        return getN0() * Math.exp(logDemographic);
    }

    private double getAnalyticIntegral(double start, double finish) {

        double integral = 0.0;
        double logDemographic = 0.0;
        double lastTransitionTime = 0.0;
        int currentEpoch = 0;

        // Account for all epochs before start
        while (currentEpoch < transitionTime.length && start > transitionTime[currentEpoch]) {
            logDemographic += -rate[currentEpoch] * (transitionTime[currentEpoch] - lastTransitionTime);
            lastTransitionTime = transitionTime[currentEpoch];
            ++currentEpoch;
        }

        // Account for all epochs before finish
        while (currentEpoch < transitionTime.length && finish > transitionTime[currentEpoch]) {

            // Add to integral
            if (rate[currentEpoch] == 0.0) {
                integral += (transitionTime[currentEpoch] - start) / Math.exp(logDemographic);
            } else {
                integral += (Math.exp(transitionTime[currentEpoch] * rate[currentEpoch]) -
                        Math.exp(start * rate[currentEpoch])) / Math.exp(logDemographic) / rate[currentEpoch];
            }

            // Update demographic function
            logDemographic += -rate[currentEpoch] * (transitionTime[currentEpoch] - lastTransitionTime);
            start = lastTransitionTime = transitionTime[currentEpoch];
            ++currentEpoch;
        }

        // End of integral
        if (rate[currentEpoch] == 0.0) {
            integral += (finish - start) / Math.exp(logDemographic);
        } else {
            integral += (Math.exp(finish * rate[currentEpoch]) -
                                    Math.exp(start * rate[currentEpoch])) / Math.exp(logDemographic) / rate[currentEpoch];
        }

        return integral / getN0();
    }

    /**
     * Calculates the integral 1/N(x) dx between start and finish.
     */
    public double getIntegral(double start, double finish) {
        // TODO An analytic solution exists
        double analytic = getAnalyticIntegral(start, finish);

        if (DEBUG) {
            double numeric = getNumericalIntegral(start, finish);

            if (Math.abs(analytic - numeric) > 1E-10) {
                System.err.println(analytic);
                System.err.println(numeric);
                throw new RuntimeException("Error in analytic calculation");
            }
        }

        return analytic;
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

    static final private boolean DEBUG = true;
}
