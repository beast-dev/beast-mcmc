/*
 * Coalescent.java
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

import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.math.Binomial;
import dr.math.MultivariateFunction;

/**
 * A likelihood function for the coalescent. Takes a tree and a demographic model.
 *
 * Parts of this class were derived from C++ code provided by Oliver Pybus.
 *
 * @version $Id: Coalescent.java,v 1.12 2005/05/24 20:25:55 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class Coalescent implements MultivariateFunction, Units {

    // PUBLIC STUFF

    public Coalescent(Tree tree, DemographicFunction demographicFunction) {
        this(new TreeIntervals(tree), demographicFunction);
    }

    public Coalescent(IntervalList intervals, DemographicFunction demographicFunction) {

        this.intervals = intervals;
        this.demographicFunction = demographicFunction;
    }


    /**
     * Calculates the log likelihood of this set of coalescent intervals,
     * given a demographic model.
     */
    public double calculateLogLikelihood() {

        return calculateLogLikelihood(intervals, demographicFunction);
    }

    /**
     * Calculates the log likelihood of this set of coalescent intervals,
     * given a demographic model.
     */
    public static double calculateLogLikelihood(IntervalList intervals, DemographicFunction demographicFunction) {
        return calculateLogLikelihood(intervals, demographicFunction, 0.0);
    }

    /**
     * Calculates the log likelihood of this set of coalescent intervals,
     * given a demographic model.
     */
    public static double calculateLogLikelihood(IntervalList intervals, DemographicFunction demographicFunction, double threshold) {

        double logL = 0.0;

        double startTime = 0.0;
        final int n = intervals.getIntervalCount();
        for (int i = 0; i < n; i++) {

            final double duration = intervals.getInterval(i);
            final double finishTime = startTime + duration;

            final double intervalArea = demographicFunction.getIntegral(startTime, finishTime);
            if( intervalArea == 0 && duration != 0 ) {
                return Double.NEGATIVE_INFINITY;
            }
            final int lineageCount = intervals.getLineageCount(i);

            final double kChoose2 = Binomial.choose2(lineageCount);
            // common part
            logL += -kChoose2 * intervalArea;

            if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {

                final double demographicAtCoalPoint = demographicFunction.getDemographic(finishTime);

                // if value at end is many orders of magnitude different than mean over interval reject the interval
                // This is protection against cases where ridiculous infinitesimal population size at the end of a
                // linear interval drive coalescent values to infinity.

                if( duration == 0.0 || demographicAtCoalPoint * (intervalArea/duration) >= threshold ) {
                    //                if( duration == 0.0 || demographicAtCoalPoint >= threshold * (duration/intervalArea) ) {
                    logL -= Math.log(demographicAtCoalPoint);
                } else {
                    // remove this at some stage
                    //  System.err.println("Warning: " + i + " " + demographicAtCoalPoint + " " + (intervalArea/duration) );
                    return Double.NEGATIVE_INFINITY;
                }

            }

            startTime = finishTime;
        }

        return logL;
    }

    /**
     * Calculates the log likelihood of this set of coalescent intervals,
     * using an analytical integration over theta.
     */
    public static double calculateAnalyticalLogLikelihood(IntervalList intervals) {

        if (!intervals.isCoalescentOnly()) {
            throw new IllegalArgumentException("Can only calculate analytical likelihood for pure coalescent intervals");
        }

        final double lambda = getLambda(intervals);
        final int n = intervals.getSampleCount();

        // assumes a 1/theta prior
        //logLikelihood = Math.log(1.0/Math.pow(lambda,n));

        // assumes a flat prior
        return (1-n) * Math.log(lambda); // Math.log(1.0/Math.pow(lambda,n-1));
    }

    /**
     * Returns a factor lambda such that the likelihood can be expressed as
     * 1/theta^(n-1) * exp(-lambda/theta). This allows theta to be integrated
     * out analytically. :-)
     */
    private static double getLambda(IntervalList intervals) {
        double lambda = 0.0;
        for (int i= 0; i < intervals.getIntervalCount(); i++) {
            lambda += (intervals.getInterval(i) * intervals.getLineageCount(i));
        }
        lambda /= 2;

        return lambda;
    }

    // **************************************************************
    // MultivariateFunction IMPLEMENTATION
    // **************************************************************

    public double evaluate(double[] argument) {
        for (int i = 0; i < argument.length; i++) {
            demographicFunction.setArgument(i, argument[i]);
        }

        return calculateLogLikelihood();
    }

    public int getNumArguments() {
        return demographicFunction.getNumArguments();
    }

    public double getLowerBound(int n) {
        return demographicFunction.getLowerBound(n);
    }

    public double getUpperBound(int n) {
        return demographicFunction.getUpperBound(n);
    }

    // **************************************************************
    // Units IMPLEMENTATION
    // **************************************************************

    /**
     * Sets the units these coalescent intervals are
     * measured in.
     */
    public final void setUnits(Type u)
    {
        demographicFunction.setUnits(u);
    }

    /**
     * Returns the units these coalescent intervals are
     * measured in.
     */
    public final Type getUnits()
    {
        return demographicFunction.getUnits();
    }

    public DemographicFunction getDemographicFunction(){
        return demographicFunction;
    }

    public IntervalList getIntervals(){
        return intervals;
    }

    /** The demographic function. */
    DemographicFunction demographicFunction = null;

    /** The intervals. */
    IntervalList intervals = null;
}