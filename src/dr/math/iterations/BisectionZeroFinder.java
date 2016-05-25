/*
 * BisectionZeroFinder.java
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
 * Zero finding by bisection.
 *
 * @author Didier H. Besset
 */
public class BisectionZeroFinder extends FunctionalIterator
{
	/**
	 * Value at which the function's value is negative.
	 */
	private double xNeg;

	/**
	 * Value at which the function's value is positive.
	 */
	private double xPos;

    /**
     * @param func DhbInterfaces.OneVariableFunction
     */
    public BisectionZeroFinder(OneVariableFunction func) {
        super(func);
    }

    /**
     * @param func DhbInterfaces.OneVariableFunction
     * @param x1 location at which the function yields a negative value
     * @param x2 location at which the function yields a positive value
     */
    public BisectionZeroFinder( OneVariableFunction func, double x1, double x2)
                                                throws IllegalArgumentException
    {
        this(func);
        setNegativeX( x1);
        setPositiveX( x2);
    }

    /**
     * @return double
     */
    public double evaluateIteration()
    {
        result = ( xPos + xNeg) * 0.5;
        if ( f.value(result) > 0 )
            xPos = result;
        else
            xNeg = result;
        return relativePrecision( Math.abs( xPos - xNeg));
    }

    /**
     * @param x double
     * @exception java.lang.IllegalArgumentException
     * 					if the function's value is not negative
     */
    public void setNegativeX( double x) throws IllegalArgumentException
    {
        if ( f.value( x) > 0 )
            throw new IllegalArgumentException( "f("+x+
                                    ") is positive instead of negative");
        xNeg = x;
    }

    /**
     * (c) Copyrights Didier BESSET, 1999, all rights reserved.
     * @param x double
     * @exception java.lang.IllegalArgumentException
     * 					if the function's value is not positive
     */
    public void setPositiveX( double x) throws IllegalArgumentException
    {
        if ( f.value( x) < 0 )
            throw new IllegalArgumentException( "f("+x+
                                    ") is negative instead of positive");
        xPos = x;
    }
}