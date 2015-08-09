/*
 * MultivariateMinimum.java
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

package dr.math;


/**
 * abstract base class for minimisation of a multivariate function
 *
 * @author Korbinian Strimmer
 */
public abstract class MultivariateMinimum
{
	//
	// Public stuff
	//

	/** total number of function evaluations necessary */
	public int numFun;

	/**
	 * maxFun is the maximum number of calls to fun allowed.
	 * the default value of 0 indicates no limit on the number
	 * of calls.
	 */
	public int maxFun = 0;

	/**
	 * numFuncStops is the number of consecutive positive
	 * evaluations of the stop criterion based on function evaluation
	 * necessary to cause the abortion of the optimization
	 * (default is 4)
	 */
	public int numFuncStops = 4;


	/**
	 * Find minimum close to vector x
	 *
	 * @param f multivariate function
	 * @param xvec initial guesses for the minimum
	 *         (contains the location of the minimum on return)
	 *
	 * @return minimal function value
	 */
	public double findMinimum(MultivariateFunction f, double[] xvec)
	{
		optimize(f, xvec, MachineAccuracy.EPSILON, MachineAccuracy.EPSILON);

		return f.evaluate(xvec);
	}
	/**
	 * Find minimum close to vector x
	 * (desired fractional digits for each parameter is specified)
	 *
	 * @param f multivariate function
	 * @param xvec initial guesses for the minimum
	 *         (contains the location of the minimum on return)
	 * @param fxFracDigits desired fractional digits in the function value
	 * @param xFracDigits desired fractional digits in parameters x
	 *
	 * @return minimal function value
	 */
	public double findMinimum(MultivariateFunction f, double[] xvec,
		int fxFracDigits, int xFracDigits)
	{
		return findMinimum(f,xvec,fxFracDigits,xFracDigits,null);
	}
	/**
	 * Find minimum close to vector x
	 * (desired fractional digits for each parameter is specified)
	 *
	 * @param f multivariate function
	 * @param xvec initial guesses for the minimum
	 *         (contains the location of the minimum on return)
	 * @param fxFracDigits desired fractional digits in the function value
	 * @param xFracDigits desired fractional digits in parameters x
	 *
	 * @return minimal function value
	 */
	public double findMinimum(MultivariateFunction f, double[] xvec,
		int fxFracDigits, int xFracDigits, MinimiserMonitor monitor)
	{
		double tolfx = Math.pow(10, -1-fxFracDigits);
		double tolx = Math.pow(10, -1-xFracDigits);

		optimize(f, xvec, tolfx, tolx,monitor);

		// trim x
		double m = Math.pow(10, xFracDigits);
		for (int i = 0;  i < xvec.length; i++)
		{
			xvec[i] = Math.round(xvec[i]*m)/m;
		}
		// trim fx
		return Math.round(f.evaluate(xvec)*m)/m;
	}

	/**
	 * The actual optimization routine
	 * (needs to be implemented in a subclass of MultivariateMinimum).
	 * It finds a minimum close to vector x when the
	 * absolute tolerance for each parameter is specified.
				 *
	 * @param f multivariate function
	 * @param xvec initial guesses for the minimum
	 *         (contains the location of the minimum on return)
	 * @param tolfx absolute tolerance of function value
	 * @param tolx absolute tolerance of each parameter
	 */
	public abstract void optimize(MultivariateFunction f, double[] xvec, double tolfx, double tolx);

	/**
	 * The actual optimization routine
	 *
	 * It finds a minimum close to vector x when the
	 * absolute tolerance for each parameter is specified.
				 *
	 * @param f multivariate function
	 * @param xvec initial guesses for the minimum
	 *         (contains the location of the minimum on return)
	 * @param tolfx absolute tolerance of function value
	 * @param tolx absolute tolerance of each parameter
	 * @param monitor A monitor object that receives information about the minimising process (for display purposes)
	 * @note The default implementation just calls the optimize function with out the Monitor!
	 */

	public void optimize(MultivariateFunction f, double[] xvec, double tolfx, double tolx, MinimiserMonitor monitor) {
		optimize(f,xvec,tolfx,tolx);
	}
	/**
	 * Checks whether optimization should stop
				 *
	 * @param fx current function value
	 * @param x current values of function parameters
	 * @param tolfx absolute tolerance of function value
	 * @param tolx absolute tolerance of each parameter
	 * @param firstCall needs to be set to true when this routine is first called
	 *        otherwise it should be set to false
	 *
	 * @return true if either x and its previous value are sufficiently similar
	 *         or if fx and its previous values are sufficiently similar
	 *         (test on function value has to be succesful numFuncStops consecutive
	 *         times)
	 */
	public boolean stopCondition(double fx, double[] x, double tolfx,
		double tolx, boolean firstCall)
	{
		boolean stop = false;

		if (firstCall)
		{
			countFuncStops = 0;
			fxold = fx;
			xold = new double[x.length];
			copy(xold, x);
		}
		else
		{
			if (xStop(x, xold, tolx))
			{
				stop = true;
			}
			else
			{
				if (fxStop(fx, fxold, tolfx))
				{
					countFuncStops++;
				}
				else
				{
					countFuncStops = 0;
				}

				if (countFuncStops >= numFuncStops)
				{
					stop = true;
				}
			}
		}

		if (!stop)
		{
			fxold = fx;
			copy(xold, x);
		}

		return stop;
	}


	/**
	 * Copy source vector into target vector
	 *
	 * @param target parameter array
	 * @param source parameter array
	 */
	public void copy(double[] target, double[] source)
	{
		for (int i = 0; i < source.length; i++)
		{
			target[i] = source[i];
		}
	}

	//
	// Private stuff
	//

	// number of fStops
	private int countFuncStops;

	// old function and parameter values
	private double fxold;
	private double[] xold;

	private boolean xStop(double[] x, double[] xold, double tolx)
	{
		boolean stop = true;

		for (int i = 0; i < x.length && stop == true; i++)
		{
			if (Math.abs(x[i]-xold[i]) > tolx)
			{
				stop = false;
			}
		}

		return stop;
	}

	private boolean fxStop(double fx, double fxold, double tolfx)
	{
		if (Math.abs(fx-fxold) > tolfx)
		{
			return false;
		}
		else
		{
			return true;
		}
	}
}
