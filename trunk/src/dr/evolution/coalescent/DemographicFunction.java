/*
 * DemographicFunction.java
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

import dr.evolution.util.Units;
import dr.math.Binomial;
import dr.math.MathUtils;

/**
 * This interface provides methods that describe a demographic function.
 *
 * Parts of this class were derived from C++ code provided by Oliver Pybus.
 *
 * @version $Id: DemographicFunction.java,v 1.12 2005/05/24 20:25:55 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Korbinian Strimmer
 */
public interface DemographicFunction extends Units {

	/**
	 * Gets the value of the demographic function N(t) at time t.
	 */
	double getDemographic(double t);

	/**
	 * Returns value of demographic intensity function at time t
	 * (= integral 1/N(x) dx from 0 to t).
	 */
	double getIntensity(double t);

	/**
	 * Returns value of inverse demographic intensity function 
	 * (returns time, needed for simulation of coalescent intervals).
	 */
	double getInverseIntensity(double x);

	/**
	 * Calculates the integral 1/N(x) dx between start and finish.
	 */
	double getIntegral(double start, double finish);
				
	/**
	 * Returns the number of arguments for this function.
	 */
	int getNumArguments();
	
	/**
	 * Returns the name of the nth argument of this function.
	 */
	String getArgumentName(int n);

	/**
	 * Returns the value of the nth argument of this function.
	 */
	double getArgument(int n);

	/**
	 * Sets the value of the nth argument of this function.
	 */
	void setArgument(int n, double value);

	/**
	 * Returns the lower bound of the nth argument of this function.
	 */
	double getLowerBound(int n);
	
	/**
	 * Returns the upper bound of the nth argument of this function.
	 */
	double getUpperBound(int n);

	/**
	 * Returns a copy of this function.
	 */
	DemographicFunction getCopy();
	
	public abstract class Abstract implements DemographicFunction
	{
			
		/**
		 * Construct demographic model with default settings
		 */
		public Abstract(int units) {
		
			setUnits(units);
		}

		// general functions

		/**
		 * Calculates the integral 1/N(x) dx between start and finish. 
		 */
		public double getIntegral(double start, double finish)
		{
			return getIntensity(finish) - getIntensity(start);
		}
			
		/**
		 * Numerically estimates the integral between start and finish.
		 */
		public double getNumericalIntegral(double start, double finish) {

			int slices = 99;
			double[] intensities = new double[slices];
			
			for (int i =0; i < intensities.length; i++) {
				double time = start + (i*finish / (slices-1));
				intensities[i] = 1.0 / getDemographic(time);
			}
			/* FIXME removed this call until we can find equivalent with compatible license
			double integral = numericalMethods.calculus.integration.NewtonCotes.simpsonSum(intensities) * (finish-start);
			return integral;
			*/
			throw new RuntimeException("Demographic model numerical integration not yet implemented");		
		}
		
	
	    // **************************************************************
	    // Units IMPLEMENTATION
	    // **************************************************************

		/**
		 * Units in which population size is measured.
		 */
		private int units;

		/**
		 * sets units of measurement.
		 *
		 * @param u units
		 */
		public void setUnits(int u)
		{
			units = u;
		}

		/**
		 * returns units of measurement.
		 */
		public int getUnits()
		{
			return units;
		}
	};

	public static class Utils
	{
		/**
		 * Returns an random interval size selected from the Kingman prior of the demographic model.
		 */
		public static double getSimulatedInterval(DemographicFunction demographicFunction, int lineageCount, double timeOfLastCoalescent)
		{
			double U = MathUtils.nextDouble(); // create unit uniform random variate
				
			double tmp = -Math.log(U)/Binomial.choose2(lineageCount) + demographicFunction.getIntensity(timeOfLastCoalescent);
			double interval = demographicFunction.getInverseIntensity(tmp) - timeOfLastCoalescent;
			
			return interval;
		}
		
		/**
		 * This function tests the consistency of the 
		 * getIntensity and getInverseIntensity methods
		 * of this demographic model. If the model is
		 * inconsistent then a RuntimeException will be thrown.
		 * @param demographicFunction the demographic model to test.
		 * @param steps the number of steps between 0.0 and maxTime to test.
		 * @param maxTime the maximum time to test.
		 */
		public static void testConsistency(DemographicFunction demographicFunction, int steps, double maxTime) {
			
			double delta = maxTime / (double)steps;
			
			for (int i =0; i <= steps; i++) {
				double time = (double)i * delta;
				double intensity = demographicFunction.getIntensity(time);
				double newTime = demographicFunction.getInverseIntensity(intensity);
							
				if (Math.abs(time-newTime) > 1e-12) {
					throw new RuntimeException(
						"Demographic model not consistent! error size = " + 
						Math.abs(time-newTime)); 
				}
			}
		}
	};
		
}