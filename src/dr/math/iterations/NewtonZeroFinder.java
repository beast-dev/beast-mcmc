/*
 * NewtonZeroFinder.java
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

package dr.math.iterations;

import dr.math.MathUtils;
import dr.math.functionEval.DrMath;
import dr.math.functionEval.FunctionDerivative;
import dr.math.interfaces.OneVariableFunction;
/**
 * Finds the zeroes of a function using Newton approximation.
 * Note: the zero of a function if the value at which the function's
 * value is zero.
 *
 * @author Didier H. Besset
 */
public class NewtonZeroFinder extends FunctionalIterator
{
	/**
	 * Derivative of the function for which the zero will be found.
	 */
	private OneVariableFunction df;
/**
* Constructor method.
* @param func the function for which the zero will be found.
* @param start the initial value for the search.
*/
public NewtonZeroFinder( OneVariableFunction func, double start)
{
	super( func);
	setStartingValue( start);
}
/**
* Constructor method.
* @param func the function for which the zero will be found.
* @param dFunc derivative of func.
* @param start the initial value for the search.
*/
public NewtonZeroFinder( OneVariableFunction func,
							OneVariableFunction dFunc, double start)
									throws IllegalArgumentException
{
	this( func, start);
	setDerivative( dFunc);
}
/**
* Evaluate the result of the current interation.
* @return the estimated precision of the result.
*/
public double evaluateIteration()
{
	double delta = f.value( result) / df.value( result);
	result -= delta;
	return relativePrecision( Math.abs( delta));
}
/**
* Initializes internal parameters to start the iterative process.
* Assigns default derivative if necessary.
*/
public void initializeIterations()
{
	if ( df == null)
		df = new FunctionDerivative( f);
	if ( Double.isNaN( result) )
		result = 0;
	int n = 0;
	while ( DrMath.equal( df.value( result), 0) )
	{
		if ( ++n > getMaximumIterations() )
			break;
		result += MathUtils.nextDouble();
	}
}
/**
 * (c) Copyrights Didier BESSET, 1999, all rights reserved.
 * @param dFunc DhbInterfaces.OneVariableFunction
 * @exception java.lang.IllegalArgumentException
 *							if the derivative is not accurate.
 */
public void setDerivative( OneVariableFunction dFunc)
										throws IllegalArgumentException
{
	df = new FunctionDerivative( f);
	if ( !DrMath.equal( df.value( result), dFunc.value( result), 0.001) )
		throw new IllegalArgumentException
						( "Supplied derative function is inaccurate");
	df = dFunc;
}
/**
 * (c) Copyrights Didier BESSET, 1999, all rights reserved.
 */
public void setFunction( OneVariableFunction func)
{
	super.setFunction( func);
	df = null;
}
/**
* Defines the initial value for the search.
*/
public void setStartingValue( double start)
{
	result = start;
}
}