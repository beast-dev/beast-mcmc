/*
 * BFloat.java
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
 *
 * BFloat -- floats with added buoyancy
 *
 * Implements numbers with precision equal to
 * floats, and range roughly 
 *
 * 1.0e-K ... 1.0e+K   with   K = log(2^(104*10^9))
 *
 * BFloat is a signed version
 * UBFloat is an unsigned version (slightly faster)
 * BDouble is a signed double version
 * UBDouble is an unsigned double version
 *
 * Todo:
 * - Implement everything except BFloat
 * - Implement comparisons
 * - Implement division, subtraction
 *
 *
 * @author Gerton Lunter
 *
 * 14/4/2003
 *
 */

 public class BFloat implements Cloneable {
 	/** Mantissa and exponent: */	
 	float iF;
 	int iE;
 	
	/** Static stuff: */	
 	/** Mantissa iF is guaranteed to satisfy iRangeInv <= |iF| <= iRange.
 	    2^104 = 2.03e+31.  Nice clean float, to minimize computational noise. */  
 	static private final float iRange = 20282409603651670423947251286016.0F;
 	static private final float iRangeInv = 1.0F/iRange;
 	static private final float iLogRange = (float)Math.log(iRange);
 	
 	/** Value between square root of iRange, and iRange.  Square of this
 	    value should still be representable with full mantissa. */
 	static private final float iHalfRange = 1.0e+18F;
 	static private final float iHalfRangeInv = 1.0e-18F;
 	
 	/** Tiniest number representable is iRangeInv ^ iInfExp */
 	static private final int iInfExp = 1000000000;
 	
 	/** true if checking for exponent overflow */
 	static private final boolean iCheckOF = false;

	/** Static array for adding numbers */
	static private final int iConvTableSize = 100;
	static private float[] iConvTable = new float[iConvTableSize];

	/** Static initialization block: */
	static {
		iConvTable[0]= 1.0F;
		for (int i=1; i<iConvTableSize; i++)
			iConvTable[i] = iConvTable[i-1] * iRangeInv;
	}


 	/** Clone this object */
    public Object clone() {
		try {
	    	// This magically creates an object of the right type
	    	BFloat iObj = (BFloat) super.clone();
	    	// No need to call clone(), as these are not Objects
	    	iObj.iE = iE;
	    	iObj.iF = iF;
	    	return iObj;
		} catch (CloneNotSupportedException e) {
	    	System.out.println("BFloat.clone: Something happened that cannot happen -- ?");
	    	return null;
		}
    }

    /** Overriding built-in toString */
    public String toString() {
		String iResult;
		if (iE == iInfExp) {
			iResult = "1.0e+Inf";
		} else if (iE == -iInfExp) {
			iResult = "1.0e-Inf";
		} else {
			double iM = (Math.log(iF) + Math.log(iRange)*(double)iE)/Math.log(10.0);
			long iExp = (long)Math.floor(iM);
			iM -= Math.floor(iM);
			float iFM = (float)Math.exp(iM*Math.log(10.0));
			// The double->float cast may result in iFM == 10.0, and we can't have that.
			if (iFM >= 10.0F) {
				iFM = 1.0F;
				iExp++;
			}
			iResult = String.valueOf(iFM) + "e";
			if (iExp<0)
				iResult += String.valueOf(iExp);
			else
				iResult += "+" + String.valueOf(iExp);
		}
		return iResult;
	}

	/** First and intermediate normalization */
 	
 	private void normalise() {
 		if (iF > 0.0) {
 			while (iF > iHalfRange) {
 				iF *= iRangeInv;
 				iE++;
 			}
 			while (iF < iHalfRangeInv) {
 				iF *= iRange;
 				iE--;
 			}
 		} else if (iF < 0.0) {
 			while (iF < -iHalfRange) {
 				iF *= iRangeInv;
 				iE++;
 			}
 			while (iF > -iHalfRangeInv) {
 				iF *= iRange;
 				iE--;
 			}
 		} else {
 			iE = -iInfExp;
 		}
 		if (iCheckOF) {
 			if (Math.abs(iE)>iInfExp) {
 				System.out.println("BFloat: exponent overflow");
 				/* throw exception */
 			}
 		}
 	}
 	
 	private void normaliseOnce() {
 		if (iF > 0.0) {
 			if (iF > iHalfRange) {
 				iF *= iRangeInv;
 				iE++;
 			} else if (iF < iHalfRangeInv) {
 				iF *= iRange;
 				iE--;
 			}
 		} else if (iF < 0.0) {
 			if (iF < -iHalfRange) {
 				iF *= iRangeInv;
 				iE++;
 			} else if (iF > -iHalfRangeInv) {
 				iF *= iRange;
 				iE--;
 			}
 		} else {
 			iE = -iInfExp;
 		}
 		if (iCheckOF) {
 			if (Math.abs(iE)>iInfExp) {
 				System.out.println("BFloat: exponent overflow");
 				/* throw exception */
 			}
 		}
 	}


	/** Various constructors */
	
 	public BFloat(float iX) {
 		iF = iX;
 		iE = 0;
 		normalise();
 	}
 	
 	public BFloat(double iX) {
 	    // This is a bit of code-duplication, but
 	    // I didn't see a clean way of doing it.
 		int iExp = 0;
  		if (iX > 0.0) {
 			while (iX > iHalfRange) {
 				iX *= iRangeInv;
 				iExp++;
 			}
 			while (iX < iHalfRangeInv) {
 				iX *= iRange;
 				iExp--;
 			}
 		} else if (iX < 0.0) {
 			while (iX < -iHalfRange) {
 				iX *= iRangeInv;
 				iExp++;
 			}
 			while (iX > -iHalfRangeInv) {
 				iX *= iRange;
 				iExp--;
 			}
 		} else {
 			iExp = -iInfExp;
 		}
 		iF = (float)iX;
 		iE = iExp;
 	}
 	
 	public BFloat(float iF, int iE) {
 		this.iF = iF;
 		this.iE = iE;
 	}

	/** Static members that return a BFloat */
	
	public static BFloat exp(float iX) {
		int iN = (int)Math.round(iX / iLogRange);
		return new BFloat( (float)Math.exp(iX - iN*iLogRange), iN );
	}
	
	/** Members that return a non-BFloat.  Nicer to overload Math.log etc, but hey. */
	
	public double log() {
		return iE*(double)iLogRange + Math.log(iF);
	}
 	
 	/** Calculations */
 	
 	public void add(BFloat iBF) {
 		if (iE >= iBF.iE) {
 			if (iE < iBF.iE + iConvTableSize) {
 				iF += iBF.iF * iConvTable[ iE - iBF.iE ];
 			}
 		} else {
 			if (iE > iBF.iE - iConvTableSize) {
 				iF = iBF.iF + iF * iConvTable[ iBF.iE - iE ];
 				iE = iBF.iE;
 			} else {
 				iF = iBF.iF;
 				iE = iBF.iE;
 			}
 		}
 		// No normalization, for speed.  Although adding numbers
 		// may break the invariant, this is would require a bit of
 		// devious programming, since iHalfRange is not precisely 
 		// the square root of iRange.
 	}
 	
 	public void add(float iX) {
 		this.add( new BFloat( iX ) );
 	}
 	
 	public void add(double iX) {
 		this.add( new BFloat( iX ) );
 	}
 	
 	public void multiply(BFloat iBF) {
 		iF *= iBF.iF;
 		iE += iBF.iE;
 		normaliseOnce();
 	}
 	
 	public void multiply(float iX) {
 		this.multiply( new BFloat( iX ) );
 	}
 	
 	public void multiply(double iX) {
 		this.multiply( new BFloat( iX ) );
 	}
 	
	/** Todo: less, lessequal,  */
	
}

