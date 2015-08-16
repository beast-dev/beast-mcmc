/*
 * DemographicFunction.java
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

import dr.evolution.util.Units;
import dr.math.Binomial;
import dr.math.MathUtils;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.MaxIterationsExceededException;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.integration.RombergIntegrator;

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
public interface DemographicFunction extends UnivariateRealFunction, Units {

    /**
     * @param t time
     * @return value of the demographic function N(t) at time t
     */
	double getDemographic(double t);

    double getLogDemographic(double t);

    /**
     * @return value of demographic intensity function at time t (= integral 1/N(x) dx from 0 to t).
     * @param t time
     */
	double getIntensity(double t);

	/**
	 * @return value of inverse demographic intensity function
	 * (returns time, needed for simulation of coalescent intervals).
	 */
	double getInverseIntensity(double x);

	/**
	 * Calculates the integral 1/N(x) dx between start and finish.
     * @param start  point
     * @param finish point
     * @return integral value
     */
	double getIntegral(double start, double finish);

	/**
	 * @return the number of arguments for this function.
	 */
	int getNumArguments();

	/**
	 * @return the name of the n'th argument of this function.
	 */
	String getArgumentName(int n);

	/**
	 * @return the value of the n'th argument of this function.
	 */
	double getArgument(int n);

	/**
	 * Sets the value of the nth argument of this function.
	 */
	void setArgument(int n, double value);

	/**
	 * @return the lower bound of the nth argument of this function.
	 */
	double getLowerBound(int n);

	/**
	 * Returns the upper bound of the nth argument of this function.
	 */
	double getUpperBound(int n);

	/**
	 * Returns a copy of this function.
	 */
//	DemographicFunction getCopy();

    /**
     * A threshold for underflow on calculation of likelihood of internode intervals.
     * Most demo functions could probably return 0.0 but (e.g.,) the Extended Skyline
     * needs a non zero value to prevent a numerical problem. 
     * @return
     */
    double getThreshold();

    public abstract class Abstract implements DemographicFunction
	{
       // private static final double LARGE_POSITIVE_NUMBER = 1.0e50;
//        private static final double LARGE_NEGATIVE_NUMBER = -1.0e50;
//        private static final double INTEGRATION_PRECISION = 1.0e-5;
//        private static final double INTEGRATION_MAX_ITERATIONS = 50;

        RombergIntegrator numericalIntegrator = null;

        /**
		 * Construct demographic model with default settings
		 */
		public Abstract(Type units) {
			setUnits(units);
        }

		// general functions

        /**
         * Default implementation
         * @param t
         * @return log(demographic at t)
         */
        public double getLogDemographic(double t) {
            return Math.log(getDemographic(t));
        }

        public double getThreshold() {
            return 0;
        }

        /**
		 * Calculates the integral 1/N(x) dx between start and finish.
		 */
		public double getIntegral(double start, double finish)
		{
			return getIntensity(finish) - getIntensity(start);
		}

        /**
         * Returns the integral of 1/N(x) between start and finish, calling either the getAnalyticalIntegral or
         * getNumericalIntegral function as appropriate.
         */
		public double getNumericalIntegral(double start, double finish) {
            // AER 19th March 2008: I switched this to use the RombergIntegrator from
            // commons-math v1.2.

            if (start > finish) {
                throw new RuntimeException("NumericalIntegration start > finish");
            }

            if (start == finish) {
                return 0.0;
            }

            if (numericalIntegrator == null) {
                numericalIntegrator = new RombergIntegrator(this);
            }

            try {
                return numericalIntegrator.integrate(start, finish);
            } catch (MaxIterationsExceededException e) {
                throw new RuntimeException(e);
            } catch (FunctionEvaluationException e) {
                throw new RuntimeException(e);
            }

//            double lastST = LARGE_NEGATIVE_NUMBER;
//            double lastS = LARGE_NEGATIVE_NUMBER;
//
//            assert(finish > start);
//
//            for (int j = 1; j <= INTEGRATION_MAX_ITERATIONS; j++) {
//                // iterate doTrapezoid() until answer obtained
//
//                double st = doTrapezoid(j, start, finish, lastST);
//                double s = (4.0 * st - lastST) / 3.0;
//
//                // If answer is within desired accuracy then return
//                if (Math.abs(s - lastS) < INTEGRATION_PRECISION * Math.abs(lastS)) {
//                    return s;
//                }
//                lastS = s;
//                lastST = st;
//            }
//
//            throw new RuntimeException("Too many iterations in getNumericalIntegral");
        }

        /**
         * Performs the trapezoid rule.
         */
//        private double doTrapezoid(int n, double low, double high, double lastS) {
//
//            double s;
//
//            if (n == 1) {
//                // On the first iteration s is reset
//                double demoLow = getDemographic(low); // Value of N(x) obtained here
//                assert(demoLow > 0.0);
//
//                double demoHigh = getDemographic(high);
//                assert(demoHigh > 0.0);
//
//                s = 0.5 * (high - low) * ( (1.0 / demoLow) + (1.0 / demoHigh) );
//            } else {
//                int it=1;
//                for (int j = 1; j < n - 1; j++) {
//                    it *= 2;
//                }
//
//                double tnm = it;	// number of points
//                double del = (high - low) / tnm;	// width of spacing between points
//
//                double x = low + 0.5 * del;
//
//                double sum = 0.0;
//                for (int j = 1; j <= it; j++) {
//                    double demoX = getDemographic(x); // Value of N(x) obtained here
//                    assert(demoX > 0.0);
//
//                    sum += (1.0 / demoX);
//                    x += del;
//                }
//                s =  0.5 * (lastS + (high - low) * sum / tnm);	// New s uses previous s value
//            }
//
//            return s;
//        }

        // **************************************************************
	    // UnivariateRealFunction IMPLEMENTATION
	    // **************************************************************

        /**
         * Return the intensity at a given time for numerical integration
         * @param x the time
         * @return the intensity
         */
        public double value(double x) {
            return 1.0 / getDemographic(x);
        }

        // **************************************************************
	    // Units IMPLEMENTATION
	    // **************************************************************

		/**
		 * Units in which population size is measured.
		 */
		private Type units;

		/**
		 * sets units of measurement.
		 *
		 * @param u units
		 */
		public void setUnits(Type u)
		{
			units = u;
		}

		/**
		 * returns units of measurement.
		 */
		public Type getUnits()
		{
			return units;
		}
	}

	public static class Utils
	{
        private static double getInterval(double U, DemographicFunction demographicFunction,
                                          int lineageCount, double timeOfLastCoalescent) {
            final double intensity = demographicFunction.getIntensity(timeOfLastCoalescent);
            final double tmp = -Math.log(U)/Binomial.choose2(lineageCount) + intensity;

            return demographicFunction.getInverseIntensity(tmp) - timeOfLastCoalescent;
        }

        private static double getInterval(double U, DemographicFunction demographicFunction, int lineageCount,
                                          double timeOfLastCoalescent, double earliestTimeOfFinalCoalescent){
            if(timeOfLastCoalescent>earliestTimeOfFinalCoalescent){
                throw new IllegalArgumentException("Given maximum height is smaller than given final coalescent time");
            }
            final double fullIntegral = demographicFunction.getIntegral(timeOfLastCoalescent,
                    earliestTimeOfFinalCoalescent);
            final double normalisation = 1-Math.exp(-Binomial.choose2(lineageCount)*fullIntegral);
            final double intensity = demographicFunction.getIntensity(timeOfLastCoalescent);

            double tmp = -Math.log(1-U*normalisation)/Binomial.choose2(lineageCount) + intensity;

            return demographicFunction.getInverseIntensity(tmp) - timeOfLastCoalescent;

        }

        /**
         * @return a random interval size selected from the Kingman prior of the demographic model.
         */
		public static double getSimulatedInterval(DemographicFunction demographicFunction,
                                                  int lineageCount, double timeOfLastCoalescent)
		{
			final double U = MathUtils.nextDouble(); // create unit uniform random variate
            return getInterval(U, demographicFunction, lineageCount, timeOfLastCoalescent);
		}

        public static double getSimulatedInterval(DemographicFunction demographicFunction, int lineageCount,
                                                  double timeOfLastCoalescent, double earliestTimeOfFirstCoalescent){
            final double U = MathUtils.nextDouble();
            return getInterval(U, demographicFunction, lineageCount, timeOfLastCoalescent,
                    earliestTimeOfFirstCoalescent);
        }

		public static double getMedianInterval(DemographicFunction demographicFunction,
                                               int lineageCount, double timeOfLastCoalescent)
		{
             return getInterval(0.5, demographicFunction, lineageCount, timeOfLastCoalescent);
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

			for (int i = 0; i <= steps; i++) {
				double time = (double)i * delta;
				double intensity = demographicFunction.getIntensity(time);
				double newTime = demographicFunction.getInverseIntensity(intensity);

				if (Math.abs(time - newTime) > 1e-12) {
					throw new RuntimeException(
						"Demographic model not consistent! error size = " +
						Math.abs(time-newTime));
				}
			}
		}
	}
}