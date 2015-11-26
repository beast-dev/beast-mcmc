/*
 * DrMath.java
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

import java.io.PrintStream;
/**
 * This class implements additional mathematical functions
 * and determines the parameters of the floating point representation.
 *
 * @author Didier H. Besset
 */
public final class DrMath
{
	/**
	 * Typical meaningful precision for numerical calculations.
	 */
	static private double defaultNumericalPrecision = 0;
	/**
	 * Typical meaningful small number for numerical calculations.
	 */
	static private double smallNumber = 0;
	/**
	 * Radix used by floating-point numbers.
	 */
	static private int radix = 0;
	/**
	 * Largest positive value which, when added to 1.0, yields 0.
	 */
	static private double machinePrecision = 0;
	/**
	 * Largest positive value which, when subtracted to 1.0, yields 0.
	 */
	static private double negativeMachinePrecision = 0;
	/**
	 * Smallest number different from zero.
	 */
	static private double smallestNumber = 0;
	/**
	 * Largest possible number
	 */
	static private double largestNumber = 0;
	/**
	 * Largest argument for the exponential
	 */
	static private double largestExponentialArgument = 0;
	/**
	 * Values used to compute human readable scales.
	 */
	private static final double scales[] = {1.25, 2, 2.5, 4, 5, 7.5, 8, 10};
	private static final double semiIntegerScales[] = {2, 2.5, 4, 5, 7.5, 8, 10};
	private static final double integerScales[] = {2, 4, 5, 8, 10};

private static void computeLargestNumber()
{
	double floatingRadix = getRadix();
	double fullMantissaNumber = 1.0d - 
						floatingRadix * getNegativeMachinePrecision();
	while ( !Double.isInfinite( fullMantissaNumber) )
	{
		largestNumber = fullMantissaNumber;
		fullMantissaNumber *= floatingRadix;
	}
}
private static void computeMachinePrecision()
{
	double floatingRadix = getRadix();
	double inverseRadix = 1.0d / floatingRadix;
	machinePrecision = 1.0d;
	double tmp = 1.0d + machinePrecision;
	while ( tmp - 1.0d != 0.0d)
	{
		machinePrecision *= inverseRadix;
		tmp = 1.0d + machinePrecision;
	}
}
private static void computeNegativeMachinePrecision()
{
	double floatingRadix = getRadix();
	double inverseRadix = 1.0d / floatingRadix;
	negativeMachinePrecision = 1.0d;
	double tmp = 1.0d - negativeMachinePrecision;
	while ( tmp - 1.0d != 0.0d)
	{
		negativeMachinePrecision *= inverseRadix;
		tmp = 1.0d - negativeMachinePrecision;
	}
}
private static void computeRadix()
{
	double a = 1.0d;
	double tmp1, tmp2;
	do { a += a;
		 tmp1 = a + 1.0d;
		 tmp2 = tmp1 - a;
		} while ( tmp2 - 1.0d != 0.0d);
	double b = 1.0d;
	while ( radix == 0)
	{
		b += b;
		tmp1 = a + b;
		radix = (int) ( tmp1 - a);
	}
}
private static void computeSmallestNumber()
{
	double floatingRadix = getRadix();
	double inverseRadix = 1.0d / floatingRadix;
	double fullMantissaNumber = 1.0d - floatingRadix * getNegativeMachinePrecision();
	while ( fullMantissaNumber != 0.0d )
	{
		smallestNumber = fullMantissaNumber;
		fullMantissaNumber *= inverseRadix;
	}
}
public static double defaultNumericalPrecision()
{
	if ( defaultNumericalPrecision == 0 )
		defaultNumericalPrecision = Math.sqrt( getMachinePrecision()); 
	return defaultNumericalPrecision;
}
/**
 * @return boolean	true if the difference between a and b is
 * less than the default numerical precision
 * @param a double
 * @param b double
 */
public static boolean equal( double a, double b)
{
	return equal( a, b, defaultNumericalPrecision());
}
/**
 * @return boolean	true if the relative difference between a and b
 * is less than precision
 * @param a double
 * @param b double
 * @param precision double
 */
public static boolean equal( double a, double b, double precision)
{
	double norm = Math.max( Math.abs(a), Math.abs( b));
	return norm < precision || Math.abs( a - b) < precision * norm;
}
public static double getLargestExponentialArgument()
{
	if ( largestExponentialArgument == 0 )
		largestExponentialArgument = Math.log(getLargestNumber());
	return largestExponentialArgument;
}
/**
 * (c) Copyrights Didier BESSET, 1999, all rights reserved.
 */
public static double getLargestNumber()
{
	if ( largestNumber == 0 )
		computeLargestNumber();
	return largestNumber;
}
public static double getMachinePrecision()
{
	if ( machinePrecision == 0 )
		computeMachinePrecision();
	return machinePrecision;
}
public static double getNegativeMachinePrecision()
{
	if ( negativeMachinePrecision == 0 )
		computeNegativeMachinePrecision();
	return negativeMachinePrecision;
}
public static int getRadix()
{
	if ( radix == 0 )
		computeRadix();
	return radix;
}
public static double getSmallestNumber()
{
	if ( smallestNumber == 0 )
		computeSmallestNumber();
	return smallestNumber;
}
public static void printParameters( PrintStream printStream)
{
	printStream.println( "Floating-point machine parameters");
	printStream.println( "---------------------------------");
	printStream.println( " ");
	printStream.println( "radix = "+ getRadix());
	printStream.println( "Machine precision = "
											+ getMachinePrecision());
	printStream.println( "Negative machine precision = "
									+ getNegativeMachinePrecision());
	printStream.println( "Smallest number = "+ getSmallestNumber());
	printStream.println( "Largest number = "+ getLargestNumber());
	return;
}
public static void reset()
{
	defaultNumericalPrecision = 0;
	smallNumber = 0;
	radix = 0;
	machinePrecision = 0;
	negativeMachinePrecision = 0;
	smallestNumber = 0;
	largestNumber = 0;
}
/**
 * This method returns the specified value rounded to
 * the nearest integer multiple of the specified scale.
 *
 * @param value number to be rounded
 * @param scale defining the rounding scale
 * @return rounded value
 */
public static double roundTo(  double value, double scale)
{
	return Math.round( value / scale) * scale;
}
	/**
	 * Round the specified value upward to the next scale value.
	 * @param the value to be rounded.
	 * @param a fag specified whether integer scale are used, otherwise double scale is used.
	 * @return a number rounded upward to the next scale value.
	 */
	public static double roundToScale( double value, boolean integerValued)
	{
		double[] scaleValues;
		int orderOfMagnitude = (int) Math.floor( Math.log( value) / Math.log( 10.0));
		if ( integerValued )
		{
			orderOfMagnitude = Math.max( 1, orderOfMagnitude);
			if ( orderOfMagnitude == 1)
				scaleValues = integerScales;
			else if ( orderOfMagnitude == 2)
				scaleValues = semiIntegerScales;
			else
				scaleValues = scales;
		}
		else
			scaleValues = scales;
		double exponent = Math.pow( 10.0, orderOfMagnitude);
		double rValue = value / exponent;
		for ( int n = 0; n < scaleValues.length; n++)
		{
			if ( rValue <= scaleValues[n])
				return scaleValues[n] * exponent;
		}
		return exponent;	// Should never reach here
	}
/**
 * (c) Copyrights Didier BESSET, 1999, all rights reserved.
 * @return double
 */
public static double smallNumber()
{
	if ( smallNumber == 0 )
		smallNumber = Math.sqrt( getSmallestNumber()); 
	return smallNumber;
}
}