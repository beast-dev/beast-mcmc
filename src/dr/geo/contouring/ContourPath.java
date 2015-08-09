/*
 * ContourPath.java
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

package dr.geo.contouring;

import dr.geo.contouring.ContourAttrib;


/**
*  <p> This object represents a single contour line or
*      path and all the data that is associated with
*      it.
*  </p>
*
*  <p>  Modified by:  Joseph A. Huwaldt  </p>
*
*  @author  Joseph A. Huwaldt   Date:  November 11, 2000
*  @version November 17, 2000
*
*
*  @author Marc Suchard
**/
public class ContourPath implements Cloneable, java.io.Serializable {

	//	Tolerance for path closure.
	private static final double kSmallX = 0.001;
	private static final double kSmallY = kSmallX;

	//	X & Y coordinate arrays.
	private double[] xArr, yArr;

	//	The level index for this contour path.
	private int levelIndex;

	//	Indicates if this path is open or closed.
	private boolean closed = false;

	//	The attributes assigned to this contour level.
	private ContourAttrib attributes;


	/**
	*  Construct a contour path or line using the given arrays
	*  of X & Y values.
	*
	*  @param  attr       Attributes assigned to this contour path.
	*  @param  levelIndex The index to then level this path belongs to.
	*  @param  x          Array of X coordinate values.
	*  @param  y          Array of Y coordinate values.
	**/
	public ContourPath(ContourAttrib attr, int levelIndex, double[] x, double[] y) {

		xArr = x;
		yArr = y;
		this.levelIndex = levelIndex;
		attributes = attr;
		int np = xArr.length;

		//	Determine if the contour path is open or closed.
		if (Math.abs(x[0] - x[np-1]) < kSmallX && Math.abs(y[0] - y[np-1]) < kSmallY) {
			closed = true;
			x[np-1] = x[0];  y[np-1] = y[0];	//	Guarantee closure.
		} else
			closed = false;						//	Contour not closed.

	}


	/**
	*  Return the X coordinate values for this contour path.
	**/
	public double[] getAllX() {
		return xArr;
	}

	/**
	*  Return the Y coordinate values for this contour path.
	**/
	public double[] getAllY() {
		return yArr;
	}

	/**
	*  Return the level index for this contour path.  The level index
	*  is an index to the level that this path belongs to:  the i'th level.
	**/
	public int getLevelIndex() {
		return levelIndex;
	}

	/**
	*  Return the attributes assigned to this contour path.
	**/
	public ContourAttrib getAttributes() {
		return attributes;
	}

	/**
	*  Returns true if this contour path is closed (loops back
	*  on itself) or false if it is not.
	**/
	public boolean isClosed() {
		return closed;
	}

	/**
	*  Make a copy of this ContourPath object.
	*
	*  @return  Returns a clone of this object.
	**/
	public Object clone() {
		ContourPath newObject = null;

		try {
			// Make a shallow copy of this object.
			newObject = (ContourPath) super.clone();

			// There is no "deep" data to be cloned.

		} catch (CloneNotSupportedException e) {
			// Can't happen.
			e.printStackTrace();
		}

		// Output the newly cloned object.
		return newObject;
	}


}
