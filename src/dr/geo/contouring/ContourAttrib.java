/*
 * ContourAttrib.java
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


/**
*  <p> This object represents the attributes assigned to a
*      contour path.  Typically, the same attributes are
*      assigned to all the contour paths of a given contour
*      level.
*  </p>
*
*  <p> Right now, the only attribute used is "level", but
*      in the future I may add more.
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
public class ContourAttrib implements Cloneable, java.io.Serializable {

	//	The level (altitude) of a contour path.
	private double level;


	/**
	*  Create a contour attribute object where only
	*  the contour level is specified.
	**/
	public ContourAttrib(double level) {
		this.level = level;
	}

	/**
	*  Return the level stored in this contour attribute.
	**/
	public double getLevel() {
		return level;
	}

	/**
	*  Set or change the level stored in this contour attribute.
	**/
	public void setLevel(double level) {
		this.level = level;
	}

	/**
	*  Make a copy of this ContourAttrib object.
	*
	*  @return  Returns a clone of this object.
	**/
	public Object clone() {
		ContourAttrib newObject = null;

		try {
			// Make a shallow copy of this object.
			newObject = (ContourAttrib) super.clone();

			// There is no "deep" data to be cloned.

		} catch (CloneNotSupportedException e) {
			// Can't happen.
			e.printStackTrace();
		}

		// Output the newly cloned object.
		return newObject;
	}

}

