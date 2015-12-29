/*
 * FunctionDerivative.java
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


import dr.math.interfaces.OneVariableFunction;
/**
 * Evaluate an approximation of the derivative of a given function.
 *
 * @author Didier H. Besset
 */
public final class FunctionDerivative implements OneVariableFunction
{
	/**
	 * Function for which the derivative is computed.
	 */
	private OneVariableFunction f;
	/**
	 * Relative interval variation to compute derivative.
	 */
	 private double relativePrecision = 0.0001;


/**
* Constructor method.
* @param func the function for which the derivative is computed.
*/
public FunctionDerivative( OneVariableFunction func)
{
	this( func, 0.000001);
}
/**
* Constructor method.
* @param func the function for which the derivative is computed.
* @param precision the relative step used to compute the derivative.
*/
public FunctionDerivative( OneVariableFunction func, double precision)
{
	f = func;
	relativePrecision = precision;
}
/**
 * Returns the value of the function's derivative
 * for the specified variable value.
 */
public double value( double x)
{
	double x1 = x == 0 ? relativePrecision
					   : x * ( 1 + relativePrecision);
	double x2 = 2 * x - x1;
	return (f.value(x1) - f.value(x2)) / (x1 - x2);
}
}