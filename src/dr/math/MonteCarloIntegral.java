/*
 * MonteCarloIntegral.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.math;

/**
 * Approximates the integral of a given function using Monte Carlo integration
 *
 * @author Alexei Drummond
 *
 */
public class MonteCarloIntegral implements Integral {

	public MonteCarloIntegral(int sampleSize) {
		this.sampleSize = sampleSize;
	}

	/**
	 * @return the approximate integral of the given function
	 * within the given range using simple monte carlo integration.
	 * @param f the function whose integral is of interest
	 * @param min the minimum value of the function
	 * @param max the  upper limit of the function
	 */
	public double integrate(UnivariateFunction f, double min, double max) {
	
		double integral = 0.0;
		
		double range = (max - min);
		for (int i =1; i <= sampleSize; i++) {
			integral += f.evaluate((MathUtils.nextDouble() * range) + min);
		}
		integral *= range/(double)sampleSize;
		return integral;
	}
	
	private int sampleSize;
}
