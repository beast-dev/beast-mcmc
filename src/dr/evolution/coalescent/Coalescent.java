/*
 * Coalescent.java
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

    public static boolean detaildPrint = false;

    /**
	 * Calculates the log likelihood of this set of coalescent intervals, 
	 * given a demographic model.
	 */
	public static double calculateLogLikelihood(IntervalList intervals, DemographicFunction demographicFunction) {
		
		double logL = 0.0;
		if( detaildPrint ) { System.err.println("new:"); }
		double startTime = 0.0;
		final int n = intervals.getIntervalCount();
		for (int i = 0; i < n; i++) {
			
			final double duration = intervals.getInterval(i);
			final double finishTime = startTime + duration;

			final double intervalArea = demographicFunction.getIntegral(startTime, finishTime);
			final int lineageCount = intervals.getLineageCount(i);

            final double kover2 = Binomial.choose2(lineageCount);
            // common part
            logL += -kover2 * intervalArea;

            if( detaildPrint ) {
                System.err.print("l = " + lineageCount + " width " + duration + " int " + intervalArea
                        + " avg " + (duration * 1/intervalArea) );


            }
            // SAMPLE or NOTHING just add likelihood for nothing to happen.
            if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {

                final double demographicAtCoalPoint = demographicFunction.getDemographic(finishTime);

                // if value at end is many orders of magnitude different than mean over interval reject this interval
                // This is protection against cases where you get  ridiculous coalescent values with infitisimal
                // population size at the end of a linear interval

                if( demographicAtCoalPoint * (intervalArea/duration) > 1e-12 ) {
                  logL += Math.log(kover2 / demographicAtCoalPoint);
                    
                  if( detaildPrint ) { System.err.println(" vatend " + demographicAtCoalPoint + " lgl " + logL); }
                } else {
                    // remove this at some stage
                    System.err.println(demographicAtCoalPoint + " " + (intervalArea/duration) );
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
	
	/** The demographic function. */
	DemographicFunction demographicFunction = null;
	
	/** The intervals. */
	IntervalList intervals = null;
}