/*
 * Trace.java
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

package dr.app.tracer;



/**
 * A simple class that stores a trace for a single statistic
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: Trace.java,v 1.11 2005/07/11 14:07:26 rambaut Exp $
 */
public class Trace {

	public static final int INITIAL_SIZE = 1000;
	public static final int INCREMENT_SIZE = 1000;
	
	public Trace(String name) {

		this.name = name;
	}
	
	public Trace(String name, double[] values) {

		this.values = new double[values.length];
		System.arraycopy(values, 0, this.values, 0, values.length);
	}
	
	/**
	 * add a value
	 */
	public void add(double value) {
	
		if (valueCount == values.length) {
			double[] newValues = new double[valueCount + INCREMENT_SIZE];
			System.arraycopy(values, 0, newValues, 0, values.length);
			values = newValues;
		}
		
		values[valueCount] = value;
		valueCount++;
	}
	
	/**
	 * add all the values in an array of doubles
	 */
	public void add(double[] values) {
		for (int i = 0; i < values.length; i++) {
			add(values[i]);
		}
	}
	
	public int getCount() { return valueCount; }
	public double getValue(int index) { return values[index]; }
	public void getValues(int start, double[] destination) { 
		System.arraycopy(values, start, destination, 0, valueCount - start);
	}
	public void getValues(int start, double[] destination, int offset) { 
		System.arraycopy(values, start, destination, offset, valueCount - start);
	}
	      
	public String getName() { return name; }
	
	//************************************************************************
	// private methods
	//************************************************************************
	
	private double[] values = new double[INITIAL_SIZE];
	private int valueCount = 0;
	private String name;
}