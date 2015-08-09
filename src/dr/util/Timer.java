/*
 * Timer.java
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

package dr.util;

public class Timer {

	private long start = 0, stop = 0;

	public void start() {
		start = System.currentTimeMillis();
	}

	public void stop() {
		stop = System.currentTimeMillis();
	}

	public void update() {
		stop = System.currentTimeMillis();
	}

	/**
	 * @return the estimated milliseconds left to complete a task.
	 */
	public long calibrate(double fraction) {
		long timeTaken = System.currentTimeMillis() - start;		
	
		return Math.round(((double)timeTaken / fraction) * (1.0 - fraction));
	}
	
	public double toSeconds() {
		update();
		return toSeconds(stop - start);
	}

	public static double toSeconds(long millis) {
		return millis / 1000.0;
	}
	
	public static double toMinutes(long millis) {
		return toSeconds(millis) / 60.0; 
	}

	public static double toHours(long millis) {
		return toMinutes(millis) / 60.0;
	}

	public static double toDays(long millis) {
		return toHours(millis) / 24.0;
	}
	
	public String toString() {
		update();
		return toString(stop - start);
	}

	public static String toString(long millis) {
		if (toDays(millis) > 1.0) return toDays(millis) + " days";
		if (toHours(millis) > 1.0) return toHours(millis) + " hours";
		if (toMinutes(millis) > 1.0) return toMinutes(millis) + " minutes";
		return toSeconds(millis) + " seconds";
	}
}
