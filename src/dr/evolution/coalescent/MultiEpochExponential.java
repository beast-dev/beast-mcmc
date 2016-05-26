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

import dr.math.matrixAlgebra.Vector;

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

    private double integrateConstant(double start, double finish, double logDemographic) {
        double integral =  (finish - start) / Math.exp(logDemographic);
        return integral;
    }

    private double integrateExponential(double start, double finish, double logDemographic, double rate) {
        double integral = (Math.exp(finish * rate) - Math.exp(start * rate)) / Math.exp(logDemographic) / rate;
//        System.err.println("\tint: " + integral + " " + start + " " + finish + " " + logDemographic + " " + rate);
//        System.err.println("\t\t" + Math.exp(finish * rate) + " - " + Math.exp(start * rate));
//        System.err.println("\t\t" + Math.exp(finish * rate - logDemographic - Math.log(rate)) + " - " + Math.exp(start * rate - logDemographic - Math.log(rate)));
        return integral;
    }

    public double getAnalyticIntegral(double start, double finish) {

        if (start == finish) {
            return 0.0;
        }

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
            double incr = 0.0;
            if (rate[currentEpoch] == 0.0) {
                integral += incr = integrateConstant(start, transitionTime[currentEpoch], logDemographic);
            } else {
                integral += incr = integrateExponential(
                        start - lastTransitionTime,
                        transitionTime[currentEpoch] - lastTransitionTime,
                        logDemographic, rate[currentEpoch]);
            }
//            System.err.println("begin incr = " + incr + " for " + start + " -> " + transitionTime[currentEpoch] + " or " +
//                    (start - lastTransitionTime) + " -> " + (transitionTime[currentEpoch] - lastTransitionTime) + " @ " + rate[currentEpoch] + " & " + Math.exp(logDemographic));

            // Update demographic function
            logDemographic += -rate[currentEpoch] * (transitionTime[currentEpoch] - lastTransitionTime);
            lastTransitionTime = transitionTime[currentEpoch];
            start = lastTransitionTime;
            ++currentEpoch;
        }

        // End of integral
        double incr = 0.0;
        if (rate[currentEpoch] == 0.0) {
            integral += incr = integrateConstant(start, finish, logDemographic);
        } else {
            integral += incr = integrateExponential(
                    start - lastTransitionTime,
                    finish - lastTransitionTime,
                    logDemographic, rate[currentEpoch]);
        }
//        System.err.println("final incr = " + incr + " for " + start + " -> " + finish + " or " +
//                (start - lastTransitionTime) + " -> " + (finish - lastTransitionTime) + " @ " + rate[currentEpoch] + " & " + Math.exp(logDemographic));


        if (Double.isNaN(integral) || Double.isInfinite(integral)) {
            System.err.println(integral + " " + start + " " + finish + new Vector(rate) + "\n");
        }

        return integral / getN0();
    }

    /**
     * Calculates the integral 1/N(x) dx between start and finish.
     */
    public double getIntegral(double start, double finish) {
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

    static final private boolean DEBUG = false;
}
