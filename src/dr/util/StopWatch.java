/*
 * StopWatch.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

public class StopWatch {

	private long start = 0, total = 0;

	public void start() {
		start = System.currentTimeMillis();
	}

	public void stop() {
		total += System.currentTimeMillis() - start;
		start = 0;
	}

	public void reset() {
		total = 0;
		start = 0;
	}

	public boolean isRunning() {
		return start > 0;
	}

	public String toString() {
		boolean running = isRunning();
		if (running) {
			stop();
		}

		String time = Timer.toString(total);

		if (running) {
			start();
		}

		return time;
	}
}
