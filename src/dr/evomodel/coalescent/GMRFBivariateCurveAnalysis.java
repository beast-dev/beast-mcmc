/*
 * GMRFBivariateCurveAnalysis.java
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

import dr.inference.trace.LogFileTraces;
import dr.inference.trace.TraceException;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

public class GMRFBivariateCurveAnalysis {

	private LogFileTraces[] lft;
	private double endTime;
	private final int[] defaultIntegrationPartitionSizes = {2,5,10,15,25,50,100,1000,5000,25000};
	private final int[] defaultDensityPartitionSizes = {10};
	
	public GMRFBivariateCurveAnalysis(String[] inputFileNames, double endTime, int burnIn) 
			throws TraceException, IOException{
		this.endTime = endTime;
		
		lft = new LogFileTraces[2];
		
		for(int i = 0; i < 2; i++){
        	File file = new File(inputFileNames[i]);
           	lft[i] = new LogFileTraces(inputFileNames[i], file);
           	lft[i].loadTraces();
           	lft[i].setBurnIn(burnIn);
        }
	}
	
	public void densityReport(){
		for(int partitionSize : defaultDensityPartitionSizes){
			densityReport(partitionSize);
		}
	}
	
	private void densityReport(int partitionSize){
		double[][] values = new double[2][];
		
		values[0] = new double[lft[0].getTraceCount()];
		values[1] = new double[lft[1].getTraceCount()];
		
		
		ArrayList<DensityResult> results = new ArrayList<DensityResult>(lft[0].getStateCount());
		
		for(int index = 0; index < lft[0].getStateCount(); index++){
			lft[0].getStateValues(index, values[0],0);
			lft[1].getStateValues(index, values[1],0);
		
			results.add(new DensityResult(getDensityReport(values, partitionSize)));
		}
				
		Collections.sort(results);
		
		double squareRootResultsSize = Math.sqrt(results.size());
		int value = (int) Math.floor(squareRootResultsSize); 
				
		double radius = results.get(value).radius;
		
		double logVolume = dr.math.MathUtils.logHyperSphereVolume(partitionSize,radius);
		
		System.out.println("Approximate Density at 0 with " + 
				partitionSize + " partitions = " + 
				1.0/(squareRootResultsSize*Math.exp(logVolume)));
		
	}
	
	private class DensityResult implements Comparable<DensityResult>{
				
		public double[] results;
		public double radius;
		
		public DensityResult(double[] results){
			this.results = results;
			
			double a = 0;
			for(double result : results){
				a += result*result;
			}
						
			radius = Math.pow(a, 1.0/results.length);
		}
		
		public int compareTo(DensityResult d){
			if(this.radius < d.radius){
				return -1;
			}else if(this.radius == d.radius){
				return 0;
			}
			
			return 1;
		}
		
		public String toString(){
			String a = "(";
			for(double result : results){
				a += Double.toString(result) + ", ";
			}
			a += " Radius = " + radius + " )\n";
			return a;
		}
		
		
	}
	
	private double[] getDensityReport(double[][] values, int partitionSize){
		double[] b = new double[partitionSize];
		
		double timeInterval = endTime/partitionSize;
		double startTime = timeInterval/2.0;
		
		for(int i = 0; i < b.length; i++){
			double[] a = {
				getPopSizeAtTime(values[0], startTime, values[0].length/2),
				getPopSizeAtTime(values[1], startTime, values[1].length/2)
			};
			startTime += timeInterval;
			
			b[i] = a[0] - a[1];
		}
		
		return b;
	}
	
	public void integrationReport(){
		System.out.println("Analyzing Difference Between Curves\n");
		for(int partitionSize : defaultIntegrationPartitionSizes){
			integrationReport(partitionSize);
		}
	}
	
	private void integrationReport(int partitionSize){
		
		double[] results = new double[lft[0].getStateCount()];
		
		double[][] values = new double[2][];
		
		values[0] = new double[lft[0].getTraceCount()];
		values[1] = new double[lft[1].getTraceCount()];
		
		double startTime = endTime/(partitionSize*2);
		double timeLength = endTime/partitionSize;
				
		for(int index = 0; index < lft[0].getStateCount(); index++){
			lft[0].getStateValues(index, values[0],0);
			lft[1].getStateValues(index, values[1],0);
		
			results[index] = getPartitionReport(values, startTime, timeLength);
		}
		
		double sum = 0;
		
		for(double a : results)
			sum += a;
		
		
		System.out.println("Partition Size " + partitionSize + " = "  + sum/results.length);
			
				
	}
	
	private double getPartitionReport(double[][] values, double startTime, double timeLength){
		double b = 0;
				
		while(startTime < endTime){
			b += getPopSizeDifference(values, startTime);
			
			startTime += timeLength;
		}
		
		b *= timeLength;
		
		return b;
	}
	
	private double getPopSizeDifference(double[][] values, double time){
		double[] a =  new double[2];
				
		for(int i = 0; i < 2; i++){
			a[i] = getPopSizeAtTime(values[i],time,values[i].length/2);
		}
		
		
		return Math.abs(a[0] - a[1]);
	}
	
	
	
	private double getPopSizeAtTime(double[] values, double time, int numberOfPopSizes){
		for(int i = 0; i < numberOfPopSizes; i++){
			if(time < values[i]){
				return values[numberOfPopSizes + i];
			}
		}
				
		return -99;
	}
	
}
