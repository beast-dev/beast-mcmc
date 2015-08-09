/*
 * GMRFPopSizeStatistic.java
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

public class GMRFPopSizeStatistic extends Statistic.Abstract{
	
	private GMRFSkyrideLikelihood gsl;
	private double[] time;

	public GMRFPopSizeStatistic(double[] time, GMRFSkyrideLikelihood gsl){
		super("Popsize");
		this.gsl = gsl;
		this.time = time;
	}
	
	public String getDimensionName(int i){
		return getStatisticName() + Double.toString(time[i]);
	}
	
	public int getDimension() {
		return time.length;
	}

	public double getStatisticValue(int dim) {
		double[] coalescentHeights = gsl.getCoalescentIntervalHeights();
		double[] popSizes = gsl.getPopSizeParameter().getParameterValues();
		
		assert popSizes.length == coalescentHeights.length;
		
		for(int i = 0; i < coalescentHeights.length; i++){
			if(time[dim] < coalescentHeights[i]){
				return popSizes[i];
			}
		}
				
		return Double.NaN;
	}


}
