/*
 * BayesianSkylinePopSizeStatistic.java
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

package dr.evomodel.coalescent;

import dr.inference.model.Statistic;

public class BayesianSkylinePopSizeStatistic extends Statistic.Abstract{

	public double time;
	public BayesianSkylineLikelihood bsl;

	public BayesianSkylinePopSizeStatistic(double time, 
			BayesianSkylineLikelihood bsl){
		this.time = time;
		this.bsl = bsl;
	}
	
	public int getDimension() {
		return 1;
	}

	public double getStatisticValue(int dim) {
		double[] heights = bsl.getGroupHeights();
		double[] sizes = bsl.getPopSizeParameter().getParameterValues();
		
		
		for(int i = 0; i < heights.length; i++){
			if(this.time < heights[i]){
				return sizes[i];
			}
		}
		
		return Double.NaN;
	}

}
