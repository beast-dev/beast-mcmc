/*
 * FunctionalIterator.java
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

import dr.math.interfaces.OneVariableFunction;

/**
 * Iterative process based on a one-variable function,
 * having a single numerical result.
 *
 * @author Didier H. Besset
 */
public abstract class FunctionalIterator extends IterativeProcess
{
	/**
	 * Best approximation of the zero.
	 */
	protected double result = Double.NaN;

    /**
	 * Function for which the zero will be found.
	 */
	protected OneVariableFunction f;

    /**
     * Generic constructor.
     * @param func OneVariableFunction
     */
    public FunctionalIterator(OneVariableFunction func)
    {
        setFunction( func);
    }

    /**
     * Returns the result (assuming convergence has been attained).
     */
    public double getResult( )
    {
        return result;
    }

    /**
     * @return double
     * @param epsilon double
     */
    public double relativePrecision( double epsilon)
    {
        return relativePrecision( epsilon, Math.abs( result));
    }

    /**
     * @param func DhbInterfaces.OneVariableFunction
     */
    public void setFunction( OneVariableFunction func)
    {
        f = func;
    }

    /**
     * @param x double
     */
    public void setInitialValue( double x)
    {
        result = x;
    }
}