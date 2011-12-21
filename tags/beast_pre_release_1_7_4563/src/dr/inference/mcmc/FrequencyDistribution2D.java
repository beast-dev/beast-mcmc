/*
 * FrequencyDistribution2D.java
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

package dr.inference.mcmc;

/**
 * This class modes a two-dimensional binned frequency distribution.
 *
 * @version $Id: FrequencyDistribution2D.java,v 1.1 2002/10/29 13:48:24 rambaut Exp $
 *
 * @author Alexei Drummond
 */
class FrequencyDistribution2D {

	double binSize1, binSize2;
	double start1 = 0.0, start2 = 0.0;
	double[][] bins;

	public FrequencyDistribution2D(int numBins, double binSize1, int numBins2, double binSize2) {
		init(numBins, binSize1, numBins2, binSize2);
	}	

	private void init(int numBins, double binSize1, int numBins2, double binSize2) {
		bins = new double[numBins][numBins2];
		this.binSize1 = binSize1; 
		this.binSize2 = binSize2; 
	}

	public void addPoint(double value1, double value2) {
	
		double diff = value1 - start1;

		int index1 = (int)(diff / binSize1);
		diff = value2 - start2;
		int index2 = (int)(diff / binSize2);
		
		if ((index1 >= 0) && (index1 < bins.length) && (index2 >= 0) && (index2 < bins[0].length)) {
			bins[index1][index2] += 1;
		}
	}

	public String toString() {
		StringBuffer buffer = new StringBuffer();
		double mean;
		buffer.append("\t");
		for (int i = 0; i < bins.length; i++) {
			mean = start2 + (binSize2 * ((double)i + 0.5));
			buffer.append(mean + "\t");
		}
		buffer.append("\n");
		for (int i = 0; i < bins.length; i++) {
			mean = start1 + (binSize1 * ((double)i + 0.5));
			buffer.append(mean + "\t");
			for (int j = 0; j < bins[i].length; j++) {
				buffer.append(bins[i][j] + "\t");
			}
			buffer.append("\n");
		}
		return new String(buffer);
	}
}
