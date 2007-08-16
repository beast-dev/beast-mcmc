/*
 * TraceDistribution.java
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

import dr.stats.DiscreteStatistics;
import dr.util.HeapSort;

/**
 * A class that stores the distribution statistics for a trace
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TraceDistribution.java,v 1.5 2005/07/11 14:07:26 rambaut Exp $
 */
public class TraceDistribution {

	public TraceDistribution(double[] values, int stepSize) {
		analyseDistribution(values, stepSize);
	}
	
	public TraceDistribution(double[] values, int stepSize, double ESS) {
		analyseDistribution(values, stepSize);
		this.ESS = ESS;
	}

	public boolean isValid() { return isValid; }

	public double getMean() { return mean; }
	public double getMedian() { return median; }
	public double getLowerHPD() { return hpdLower; }
	public double getUpperHPD() { return hpdUpper; }
	public double getLowerCPD() { return cpdLower; }
	public double getUpperCPD() { return cpdUpper; }
	public double getESS() { return ESS; }
	
	public double getMinimum() { return minimum; }
	public double getMaximum() { return maximum; }
	/**
	 * Analyze trace
	 */
	private void analyseDistribution(double[] values, int stepSize) {
	
		mean = DiscreteStatistics.mean(values);

		minimum = Double.POSITIVE_INFINITY;
		maximum = Double.NEGATIVE_INFINITY;
		
		for (int i = 0; i < values.length; i++) {
			if (values[i] < minimum) minimum = values[i];
			if (values[i] > maximum) maximum = values[i];
		}
		
		if (maximum == minimum) {
			isValid = false;
			return;
		}
		
		int[] indices = new int[values.length];
		HeapSort.sort(values, indices);
		median = DiscreteStatistics.quantile(0.5, values, indices);
		cpdLower = DiscreteStatistics.quantile(0.025, values, indices);
		cpdUpper = DiscreteStatistics.quantile(0.975, values, indices);
		calculateHPDInterval(0.95, values, indices);
		ESS = values.length;
		
		isValid = true;
	}
		
	/**
	 * @param proportion the proportion of probability mass oncluded within interval.
	 */
	private void calculateHPDInterval(double proportion, double[] array, int[] indices) {
	
		double minRange = Double.MAX_VALUE;
		int hpdIndex = 0;
	
		int diff = (int)Math.round(proportion * (double)array.length);
		for (int i =0; i <= (array.length - diff); i++) {
			double minValue = array[indices[i]];
			double maxValue = array[indices[i+diff-1]];
			double range = Math.abs(maxValue - minValue);
			if (range < minRange) {
				minRange = range;
				hpdIndex = i;
			}
		}
		hpdLower = array[indices[hpdIndex]];
		hpdUpper = array[indices[hpdIndex+diff-1]];
	}

	//************************************************************************
	// private methods
	//************************************************************************
	
	protected boolean isValid = false;
	
	protected double minimum, maximum;
	protected double mean, median;
	protected double cpdLower, cpdUpper, hpdLower, hpdUpper;
	protected double ESS;
}