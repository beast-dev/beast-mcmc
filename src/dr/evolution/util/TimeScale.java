/*
 * TimeScale.java
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

package dr.evolution.util;

/**
 * Class for defining and converting between time scales
 *
 * @version $Id: TimeScale.java,v 1.11 2005/05/24 20:25:57 rambaut Exp $
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class TimeScale implements Units { 

	/**
	 * Constructor for a timescale that starts at Jan 1st 1970. This constructor
	 * can be used if we are not interested in absolute time.
	 * param units The units
	 * param backwards True if the timescale goes backwards in time
	 */
	public TimeScale(Type units, boolean backwards) {
		this(units, backwards, 0.0);
	}
	
	/**
	 * Constructor for a timescale that starts at origin
	 * param units The units
	 * param backwards True if the timescale goes backwards in time
	 * param origin The origin specified relative to 1970 in same units
	 */
	public TimeScale(Type units, boolean backwards, double origin) {
		this.units = units;
		this.backwards = backwards;
		this.origin = origin;
	}
	
	/**
	 * Constructor for a timescale that starts at origin
	 * param units The units
	 * param backwards True if the timescale goes backwards in time
	 * param origin The origin specified as a date
	 */
	public TimeScale(Type units, boolean backwards, java.util.Date origin) {
		this.units = units;
		this.backwards = backwards;
		
		long millisAhead = origin.getTime();
		
		double daysAhead = ((double)millisAhead)/MILLIS_PER_DAY;
		
		switch (units) {	
			case DAYS: this.origin = daysAhead; break;
			case MONTHS: this.origin = daysAhead / DAYS_PER_MONTH; break;
			case YEARS: this.origin = daysAhead / DAYS_PER_YEAR; break;
			default: throw new IllegalArgumentException();
		}
		
	}
	
	/**
	 * @return the units of this timescale
	 */
	public Type getUnits() { return units; }
	
	/**
	 * Sets the units for this timescale.
	 */
	public void setUnits(Type units) { this.units = units; }
	
	/**
	 * @return true if larger numbers represent older dates in this time scale.
	 */
	public boolean isBackwards() { return backwards; }
	
	/**
	 * @return the date corresponding to the value zero in this time scale.
	 */
	public double getOrigin() { return origin; }
	
	/** 
	 * @return a time given in timescale scale as a time in this timescale
	 */
	public double convertTime(double time, TimeScale timeScale) {
		
		// make it forwards
		if (timeScale.isBackwards()) time = -time;
		
		// make it absolute
		time += timeScale.getOrigin();
		
		// convert to the new timescale units
		double newTime = convertTimeUnits(time, getUnits(), timeScale.getUnits());
		
		// make it relative
		newTime -= origin;
		
		// make it backwards if required
		if (backwards) newTime = -newTime;

		return newTime;
	}
	
	public String toString() {
	
		StringBuffer buffer = new StringBuffer("timescale(");
		buffer.append(unitString(0.0));
		if (backwards) {
			buffer.append(", backwards");
		} else {
			buffer.append(", forewards");
		}
		buffer.append(" from " + origin + ")");
		
		return buffer.toString();
	}
	
	public String unitString(double time) {
		String unitString = null;
		switch (units) {	
			case DAYS: unitString = "day"; break;
			case MONTHS: unitString = "month"; break;
			case YEARS: unitString = "year"; break;
			default: throw new IllegalArgumentException();
		}
		if (time == 1.0) {
			return unitString;
		} else return unitString + "s";
	}
	
	public static void main(String[] args) {
	
		TimeScale timeScale1 = new TimeScale(Units.Type.DAYS, true);
		TimeScale timeScale2 = new TimeScale(Units.Type.YEARS, true);
		
		System.out.println(timeScale1);
		System.out.println(timeScale2);
		
		double testTime = 100.0;
		System.out.println("Test time = " + testTime);
		
		System.out.println("timeScale1.convertTime(" + testTime + ", timeScale2)=" + timeScale1.convertTime(testTime, timeScale2));
		System.out.println("timeScale2.convertTime(" + testTime + ", timeScale1)=" + timeScale2.convertTime(testTime, timeScale1));
	}
	
	
	//*************************************************************************
	// STATIC STUFF
	//*************************************************************************
	
	/** 
	 * @return time in currentUnits as newUnits.
	 */
	public static double convertTimeUnits(double time, Type currentUnits, Type newUnits) {
		
		return time * getScale(currentUnits, newUnits);
	}
	
	/** 
	 * @return the scaling factor for converting currentUnits into newUnits.
	 */
	public static double getScale(Type currentUnits, Type newUnits) {
		if (currentUnits == newUnits) return 1.0;
		
		switch (currentUnits) {
			case DAYS:
				switch (newUnits) {
					case MONTHS: return 1.0/DAYS_PER_MONTH;
					case YEARS: return 1.0/DAYS_PER_YEAR;
					default: throw new IllegalArgumentException();
				}
			case MONTHS:
				switch (newUnits) {
					case DAYS: return DAYS_PER_MONTH;
					case YEARS: return 1.0/MONTHS_PER_YEAR;
					default: throw new IllegalArgumentException();
				}
			case YEARS:
				switch (newUnits) {
					case DAYS: return DAYS_PER_YEAR;
					case MONTHS: return MONTHS_PER_YEAR;
					default: throw new IllegalArgumentException();
				}
			default: throw new IllegalArgumentException();
		}
	}
	
	//*************************************************************************
	// PRIVATE STUFF
	//*************************************************************************

	// The origin is specified in days relative to 1st January 1970
	protected double origin = 720035.0;
	protected Type units;
	protected boolean backwards;
	
	protected static double MILLIS_PER_DAY = 86400000.0;
	protected static double DAYS_PER_YEAR = 365.25;
	protected static double MONTHS_PER_YEAR = 12.0;
	protected static double DAYS_PER_MONTH = DAYS_PER_YEAR / MONTHS_PER_YEAR;
}
