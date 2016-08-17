/*
 * IntMathVec.java
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

package dr.oldevomodel.indel;

/**
 *
 *  IntMathVec -- mathematical vector of integers.
 *
 *  Vectors provide equals() and hashCode() based on value-identity,
 *  so can be used in maps, but are NOT immutable, so take care not
 *  to change the value of map-keys.
 *
 *  @author Gerton Lunter
 *  
 *  20/3/2003
 *
 */

class IntMathVec implements Cloneable {
    public int[] iV;

    /** Constructor, null vector of given length */
    IntMathVec(int iLen) {
	iV = new int[iLen];
    }

    /** Constructor, copy of existing vector */
    IntMathVec(IntMathVec iVec) {
	iV = iVec.iV.clone();
    }

    /** Constructor, copy of integer array - handy for printing */
    IntMathVec(int[] iArr) {
	iV = iArr.clone();
    }

    /** Internal, to make sure vectors are of same length. */
    // I know, I should throw an exception..  Will likely happen automatically.
    private void check(IntMathVec iVec) {
	if (iVec.iV.length != iV.length) {
	    System.out.println("IntMathVec.check: Vector sizes don't match.");
	}
    }

    /** Make equals() reflect value-identity instead of object-identity */
    public boolean equals(Object iObj) {
	if (iObj instanceof IntMathVec) {
	    IntMathVec iVec = (IntMathVec)iObj;
	    if (iVec.iV.length != iV.length)
		return false;
	    for (int i=0; i<iV.length; i++) {
		if (iV[i] != iVec.iV[i])
		    return false;
	    }
	    return true;
	}
	return false;
    }

    /** hashCode() reflects value-identity instead of object-identity */
    public int hashCode() {
        int iCode = 0;
        for( int anIV : iV ) {
            iCode = (iCode * 75) + anIV;
        }
        return iCode;
    }

    /** Clone this object */
    public IntMathVec clone() {
        try {
            // This magically creates an object of the right type
            IntMathVec iObj = (IntMathVec) super.clone();
            iObj.iV = iV.clone();
            return iObj;
        } catch (CloneNotSupportedException e) {
            System.out.println("IntMathVec.clone: Something happened that cannot happen -- ?");
            return null;
        }
    }

    /** Overriding built-in toString, produces Mathematica-readable output */
    public String toString() {
	String iResult = "{";
	for (int i=0; i<iV.length; i++) {
	    if (i != 0)
		iResult += ",";
	    iResult += iV[i];
	}
	iResult += "}";
	return iResult;
    }

    /** Now the math thingies */
    public int innerProduct(IntMathVec iVec) {
	check(iVec);
	int iSum = 0;
	for (int i=0; i<iV.length; i++)
	    iSum += iVec.iV[i]*iV[i];
	return iSum;
    }

    public boolean zeroEntry() {
        for(int anIV : iV) {
            if( anIV == 0 ) {
                return true;
            }
        }
	return false;
    }

    public void assign(IntMathVec iVec) {
	iV = iVec.iV.clone();
    }

    public void add(IntMathVec iVec) {
	check(iVec);
	for (int i=0; i<iV.length; i++)
	    iV[i] += iVec.iV[i];
    }

    public void addMultiple(IntMathVec iVec, int iMultiple) {
	check(iVec);
	for (int i=0; i<iV.length; i++)
	    iV[i] += iVec.iV[i] * iMultiple;
    }

    public void subtract(IntMathVec iVec) {
	check(iVec);
	for (int i=0; i<iV.length; i++)
	    iV[i] -= iVec.iV[i];
    }
    
}
	


