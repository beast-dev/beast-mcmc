/*
 * FrequencyDistribution.java
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

/**
 * @author Alexei Drummond
 * @version $Id: FrequencyDistribution.java,v 1.3 2005/05/24 20:26:01 rambaut Exp $
 */
public class FrequencyDistribution {

	/** The size of the bin */
	private double binSize;
	
	/** the lower boundary of first bin */
	private double start = 0.0;
	
	/** the number of values lower than first bin */
	private double smaller;

	/** the number of values greater than last bin */
	private double larger;

	/** the frequncy counts of each bin */
	private int[] bins;

	public FrequencyDistribution(double start, int numBins, double binSize) {
		init(start, numBins, binSize);
	}

	public FrequencyDistribution(int numBins, double binSize) {
		init(0.0, numBins, binSize);
	}

	public FrequencyDistribution(double[] stats, int numBins, double binSize) {
		init(0.0, numBins, binSize);
        for (double stat : stats) {
            addValue(stat);
        }
    }

	/**
	 * Returns the number of bins.
	 */
	public int getBinCount() {
		return bins.length;
	}

	public double getBinSize() {
		return binSize;
	}

	/** Returns lower bound of first bin. */
	public double getLowerBound() {
		return start;
	}

	/**
	 * @return the number of values falling in this bin.
	 */
	public int getFrequency(int bin) {
		return bins[bin];
	}

    public double getProbability(int bin) {
        int total = 0;
        for (int b : bins) {
            total = total + b;
        }
        if (total == 0) {
            return 0.0;
        } else {
		    return (double) bins[bin] / (double) total;
        }
	}

	public void addValue(double value) {
	
		double diff = value - start;

		int index = (int)(diff / binSize);
		if (index < 0) {
			smaller += 1;
		} else if (index >= bins.length) {
			larger += 1;
		} else bins[index] += 1;
	}

	public String toString() {
		StringBuffer buffer = new StringBuffer();
        buffer.append("< ").append(start).append("\t").append(smaller).append("\n");
		double mean;
		for (int i = 0; i < bins.length; i++) {
			mean = start + (binSize * ((double)i + 0.5));
            buffer.append(mean).append("\t").append(bins[i]).append("\n");
		}
		double end = start + (binSize * bins.length);
        buffer.append(">= ").append(end).append("\t").append(larger).append("\n");
		return new String(buffer);
	}

	private void init(double start, int numBins, double binSize) {
		bins = new int[numBins];
		this.binSize = binSize; 
		smaller = 0;
		larger = 0;
		this.start = start;
	}
}
