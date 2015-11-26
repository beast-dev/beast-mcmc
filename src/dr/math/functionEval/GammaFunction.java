/*
 * GammaFunction.java
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

package dr.math.functionEval;

/**
 * Gamma function (Euler's integral).
 * 
 * @author Didier H. Besset
 */
public final class GammaFunction
{
	static double sqrt2Pi = Math.sqrt( 2 * Math.PI);
	static double[] coefficients = { 76.18009172947146,
								    -86.50532032941677,
								     24.01409824083091,
								     -1.231739572450155,
								      0.1208650973866179e-2,
								     -0.5395239384953e-5};

/**
 * @return double		beta function of the arguments
 * @param x double
 * @param y double
 */
public static double beta ( double x, double y)
{
	return Math.exp( logGamma( x) + logGamma( y) - logGamma( x + y));
}
/**
 * @return long	factorial of n
 * @param n long
 */
public static long factorial ( long n)
{
	return n < 2 ? 1 : n * factorial( n - 1);
}
/**
 * @return double		gamma function
 * @param x double
 */
public static double gamma ( double x)
{
	return x > 1
				? Math.exp( leadingFactor(x)) * series(x) * sqrt2Pi / x
				: ( x > 0 ? gamma(x + 1) / x
								: Double.NaN);
}
/**
 * @return double
 * @param x double
 */
private static double leadingFactor ( double x)
{
	double temp = x + 5.5;
	return Math.log( temp) * ( x + 0.5) - temp;
}
/**
 * @return double	logarithm of the beta function of the arguments
 * @param x double
 * @param y double
 */
public static double logBeta ( double x, double y)
{
	return logGamma( x) + logGamma( y) - logGamma( x + y);
}
/**
 * @return double		log of the gamma function
 * @param x double
 */
public static double logGamma ( double x)
{
	return x > 1
				? leadingFactor(x) + Math.log( series(x) * sqrt2Pi / x)
				: ( x > 0 ? logGamma(x + 1) - Math.log( x)
								: Double.NaN);
}
/**
 * @return double		value of the series in Lanczos formula.
 * @param x double
 */
private static double series( double x)
{
	double answer = 1.000000000190015;
	double term = x;
	for ( int i = 0; i < 6; i++)
	{
		term += 1;
		answer += coefficients[i] / term;
	}	
	return answer;
}
}