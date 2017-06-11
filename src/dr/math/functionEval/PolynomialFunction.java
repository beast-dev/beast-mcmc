/*
 * PolynomialFunction.java
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

import dr.math.MathUtils;
import dr.math.interfaces.OneVariableFunction;
import dr.math.iterations.NewtonZeroFinder;

import java.util.Enumeration;
import java.util.Vector;

/**
 * Mathematical polynomial:
 * c[0] + c[1] * x + c[2] * x^2 + ....
 *
 * @author Didier H. Besset
 */
public class PolynomialFunction implements OneVariableFunction
{
	/**
	 * Polynomial coefficients.
	 */
	private final double[] coefficients;


/**
 * Constructor method.
 * @param coeffs polynomial coefficients.
 */
public PolynomialFunction( double[] coeffs)
{
	coefficients = coeffs;
}
/**
 *
 * @param r double	number added to the polynomial.
 * @return DhbFunctionEvaluation.PolynomialFunction
 */
public PolynomialFunction add( double r)
{
	int n = coefficients.length;
	double coef[] = new double[n];
	coef[0] = coefficients[0] + r;
	for ( int i = 1; i < n; i++)
		coef[i] = coefficients[i];
	return new PolynomialFunction( coef);
}
/**
 *
 * @param p DhbFunctionEvaluation.PolynomialFunction
 * @return DhbFunctionEvaluation.PolynomialFunction
 */
public PolynomialFunction add( PolynomialFunction p)
{
	int n = Math.max( p.degree(), degree()) + 1;
	double[] coef = new double[n];
	for ( int i = 0; i < n; i++ )
		coef[i] = coefficient(i) + p.coefficient(i);
	return new PolynomialFunction( coef);
}
/**
 * Returns the coefficient value at the desired position
 * @param n int	 the position of the coefficient to be returned
 * @return double the coefficient value
 */
public double coefficient( int n)
{
	return n < coefficients.length ? coefficients[n] : 0;
}
/**
 * 
 * @param r double	 a root of the polynomial (no check made).
 * @return PolynomialFunction the receiver divided by polynomial (x - r).
 */
public PolynomialFunction deflate( double r)
{
	int n = degree();
	double remainder = coefficients[n];
	double[] coef = new double[n];
	for ( int k = n - 1; k >= 0; k--)
	{
		coef[k] = remainder;
		remainder = remainder * r + coefficients[k];
	}
	return new PolynomialFunction( coef);
}
/**
 * Returns degree of this polynomial function
 * @return int degree of this polynomial function
 */
public int degree()
{
	return coefficients.length - 1;
}
/**
 * Returns the derivative of the receiver.
 * @return PolynomialFunction derivative of the receiver.
 */
public PolynomialFunction derivative()
{
	int n = degree();
	if ( n == 0 )
	{
		double coef[] = {0};
		return new PolynomialFunction( coef);
	}
	double coef[] = new double[n];
	for ( int i = 1; i <= n; i++)
		coef[i-1] = coefficients[i]*i;
	return new PolynomialFunction( coef);
}
/**
 *
 * @param r double
 * @return DhbFunctionEvaluation.PolynomialFunction
 */
public PolynomialFunction divide( double r)
{
	return multiply( 1 / r);
}
/**
 *
 * @param p PolynomialFunction
 * @return PolynomialFunction
 */
public PolynomialFunction divide( PolynomialFunction p)
{
	return divideWithRemainder(p)[0];
}
/**
 *
 * @param p PolynomialFunction
 * @return DhbFunctionEvaluation.PolynomialFunction
 */
public PolynomialFunction[] divideWithRemainder( PolynomialFunction p)
{
	PolynomialFunction[] answer = new PolynomialFunction[2];
	int m = degree();
	int n = p.degree();
	if ( m < n )
	{
		double[] q = {0};
		answer[0] = new PolynomialFunction( q);
		answer[1] = p;
		return answer;
	}
	double[] quotient = new double[ m - n + 1];
	double[] coef = new double[ m + 1];
	for ( int k = 0; k <= m; k++ )
		coef[k] = coefficients[k];
	double norm = 1 / p.coefficient( n);
	for ( int k = m - n; k >= 0; k--)
	{
		quotient[k] = coef[ n + k] * norm;
		for ( int j = n + k - 1; j >= k; j--)
			coef[j] -= quotient[k] * p.coefficient(j-k);
	}
	double[] remainder = new double[n];
	for ( int k = 0; k < n; k++)
		remainder[k] = coef[k];
	answer[0] = new PolynomialFunction( quotient);
	answer[1] = new PolynomialFunction( remainder);
	return answer;
}
/**
 * Returns the integral of the receiver having the value 0 for X = 0.
 * @return PolynomialFunction integral of the receiver.
 */
public PolynomialFunction integral( )
{
	return integral( 0);
}
/**
 * Returns the integral of the receiver having the specified value for X = 0.
 * @param value double	value of the integral at x=0
 * @return PolynomialFunction integral of the receiver.
 */
public PolynomialFunction integral( double value)
{
	int n = coefficients.length + 1;
	double coef[] = new double[n];
	coef[0] = value;
	for ( int i = 1; i < n; i++)
		coef[i] = coefficients[i-1]/i;
	return new PolynomialFunction( coef);
}
/**
 *
 * @param r double
 * @return DhbFunctionEvaluation.PolynomialFunction
 */
public PolynomialFunction multiply( double r)
{
	int n = coefficients.length;
	double coef[] = new double[n];
	for ( int i = 0; i < n; i++)
		coef[i] = coefficients[i] * r;
	return new PolynomialFunction( coef);
}
/**
 *
 * @param p DhbFunctionEvaluation.PolynomialFunction
 * @return DhbFunctionEvaluation.PolynomialFunction
 */
public PolynomialFunction multiply( PolynomialFunction p)
{
	int n = p.degree() + degree();
	double[] coef = new double[n + 1];
	for ( int i = 0; i <= n; i++)
	{
		coef[i] = 0;
		for ( int k = 0; k <= i; k++)
			coef[i] += p.coefficient(k) * coefficient(i-k);
	}
	return new PolynomialFunction( coef);
}
/**
 * 
 * @return double[]
 */
public double[] roots()
{
	return roots( DrMath.defaultNumericalPrecision());
}
/**
 *
 * @param desiredPrecision double
 * @return double[]
 */
public double[] roots( double desiredPrecision)
{
	PolynomialFunction dp = derivative();
	double start = 0;
	while ( Math.abs( dp.value( start)) < desiredPrecision )
		start = MathUtils.nextDouble();
	PolynomialFunction p = this;
	NewtonZeroFinder rootFinder = new NewtonZeroFinder( this, dp, start);
	rootFinder.setDesiredPrecision( desiredPrecision);
	Vector rootCollection = new Vector( degree());
	while ( true)
	{
		rootFinder.evaluate();
		if ( !rootFinder.hasConverged() )
			break;
		double r = rootFinder.getResult();
		rootCollection.addElement(r);
		p = p.deflate( r);
		if ( p.degree() == 0 )
			break;
		rootFinder.setFunction( p);
		try { rootFinder.setDerivative( p.derivative());}
				catch ( IllegalArgumentException e) {}
    }
	double[] roots = new double[ rootCollection.size()];
	Enumeration e = rootCollection.elements();
	int n = 0;
	while ( e.hasMoreElements() )
	{
		roots[n++] = (Double) e.nextElement();
	}
	return roots;
}
/**
 *
 * @param r double
 * @return PolynomialFunction
 */
public PolynomialFunction subtract( double r)
{
	return add( -r);
}
/**
 * 
 * @return PolynomialFunction
 * @param p PolynomialFunction
 */
public PolynomialFunction subtract( PolynomialFunction p)
{
	int n = Math.max( p.degree(), degree()) + 1;
	double[] coef = new double[n];
	for ( int i = 0; i < n; i++ )
		coef[i] = coefficient(i) - p.coefficient(i);
	return new PolynomialFunction( coef);
}
/**
* Returns a string representing the receiver
*/
public String toString()
{
	StringBuffer sb = new StringBuffer();
	boolean firstNonZeroCoefficientPrinted = false;
	for ( int n = 0; n < coefficients.length; n++)
	{
		if ( coefficients[n] != 0 )
		{
			if ( firstNonZeroCoefficientPrinted)
				sb.append( coefficients[n] > 0 ? " + " : " ");
			else
				firstNonZeroCoefficientPrinted = true;
			if ( n == 0 || coefficients[n] != 1)
				sb.append( Double.toString( coefficients[n]) );
			if ( n > 0 )
				sb.append( " X^"+n);
		}
	}
	return sb.toString();
}
/**
 * Returns the value of the polynomial for the specified variable value.
 * @param x double	value at which the polynomial is evaluated
 * @return double polynomial value.
 */
public double value( double x)
{
	int n = coefficients.length;
	double answer = coefficients[--n];
	while ( n > 0 )
		answer = answer * x + coefficients[--n];
	return answer;
}
/**
 * Returns the value and the derivative of the polynomial 
 * for the specified variable value in an array of two elements
 * @param x double	value at which the polynomial is evaluated
 * @return double[0]   the value of the polynomial
 * @return double[1]   the derivative of the polynomial
 */
public double[] valueAndDerivative( double x)
{
	int n = coefficients.length;
	double[] answer = new double[2];
	answer[0] = coefficients[--n];
	answer[1] = 0;
	while ( n > 0 )
	{
		answer[1] = answer[1] * x + answer[0];
		answer[0] = answer[0] * x + coefficients[--n];
	}
	return answer;
}
}