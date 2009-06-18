/*
 * IntersectionBounds.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.inference.model;

import java.util.ArrayList;


public class IntersectionBounds implements Bounds {
	
	IntersectionBounds(int dimension) {
		this.dimension = dimension;
	}
	
	public void addBounds(Bounds boundary) {
		if (boundary.getBoundsDimension() != dimension) {
			throw new IllegalArgumentException("Incorrect dimension of bounds, expected " + 
				dimension + " but received " + boundary.getBoundsDimension());
		}
		if (bounds == null) { 
			bounds = new ArrayList<Bounds>();
		}
		bounds.add(boundary);
	}
		
	/**
	 * Gets the maximum lower limit of this parameter and all its slave parameters
	 */
	public double getLowerLimit(int index) { 
		
		double lower = Double.NEGATIVE_INFINITY; 
		if (bounds != null) {
            for (Bounds boundary : bounds) {
                if (boundary.getLowerLimit(index) > lower) {
                    lower = boundary.getLowerLimit(index);
                }
            }
        }
		return lower;
	}
			
	/**
	 * Gets the minimum upper limit of this parameter and all its slave parameters
	 */
	public double getUpperLimit(int index) { 
		
		double upper = Double.POSITIVE_INFINITY;
		if (bounds != null) {
            for (Bounds boundary : bounds) {
                if (boundary.getUpperLimit(index) < upper) {
                    upper = boundary.getUpperLimit(index);
                }
            }
        }
		return upper;
	}
		
	public int getBoundsDimension() { return dimension; }
	
	public String toString() {
			String str = "upper=["+getUpperLimit(0);
			for (int i = 1; i < getBoundsDimension(); i++) {
				str += ", " + getUpperLimit(i);
			}
			str += "] lower=["+getLowerLimit(0);
			for (int i = 1; i < getBoundsDimension(); i++) {
				str += ", " + getLowerLimit(i);
			}
			
			str += "]";
			return str;
		}
		
	private ArrayList<Bounds> bounds = null;
	private int dimension;
}
