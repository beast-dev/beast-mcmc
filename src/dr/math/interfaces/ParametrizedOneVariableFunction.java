/*
 * ParametrizedOneVariableFunction.java
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

package dr.math.interfaces;

/**
 * ParametrizedOneVariableFunction is an interface for mathematical
 * functions of one variable depending on several parameters,
 * that is functions of the form f(x;p), where p is a vector.
 *
 * @author Didier H. Besset
 */
public interface ParametrizedOneVariableFunction
										extends OneVariableFunction
{
/**
 * @return double[]	array containing the parameters
 */
double[] parameters();
/**
 * @param p double[]	assigns the parameters
 */
void setParameters( double[] p);
/**
 * Evaluate the function and the gradient of the function with respect
 * to the parameters.
 * @return double[]	0: function's value, 1,2,...,n function's gradient
 * @param x double
 */
double[] valueAndGradient( double x);
}