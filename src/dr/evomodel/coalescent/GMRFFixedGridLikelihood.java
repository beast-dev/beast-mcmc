/*
 * GMRFFixedGridLikelihood.java
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

import java.util.ArrayList;
import java.util.logging.Logger;

import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.math.Binomial;

public class GMRFFixedGridLikelihood extends GMRFSkyrideLikelihood{

	private Parameter covariateData;
	private Parameter covariateTimes;

	private ArrayList<CoalescentIntervalWithData> intervals;
	private ArrayList<CoalescentIntervalWithData> storedIntervals;

	public static void main(String[] args){
		
		try{
			run();
		}catch(Exception e){
			System.err.println(e.getMessage());
		}
	}
	
	public static void run() throws Exception{
		NewickImporter importer = new NewickImporter("((((5:0.5,1:0.2):0.5,0:1):0.2,2:0.8):0.2,3:1.4)");
		Tree tree = importer.importNextTree();
		
		double[] data = new double[15];
		double[] times = new double[15];
		
		data[0] = 1.0;
		times[0] = 0.05;
		
		for(int i = 1; i < data.length; i++){
			data[i] = data[i-1] + 0.5;
			times[i] = times[i-1] + 0.1;
		}
		
		GMRFFixedGridLikelihood like = new GMRFFixedGridLikelihood(tree,
				new Parameter.Default(data),new Parameter.Default(times),4);
		
		System.out.println(like.getLogLikelihood());
	}
	
	public GMRFFixedGridLikelihood(Tree tree, Parameter data, Parameter times, int tips){
		super(tree, new Parameter.Default(tips), null, new Parameter.Default(5.0),
				new Parameter.Default(1.0), null, null,false, true);
				
		covariateData = data;
		covariateTimes = times;
		
		fieldLength += covariateData.getDimension();
		
		intervals = new ArrayList<CoalescentIntervalWithData>(fieldLength);
		storedIntervals = new ArrayList<CoalescentIntervalWithData>(fieldLength);
		
		sSetupIntervals();
	}
	
	public void initializationReport() {
		
	}

	public GMRFFixedGridLikelihood(Tree tree, Parameter popParameter, Parameter precParameter,
	                                      Parameter lambda, Parameter beta, MatrixParameter dMatrix,
	                                      Parameter data, Parameter times) {
		super(tree, popParameter, null, precParameter, lambda, beta, dMatrix, false, true);

		covariateData = data;
		covariateTimes = times;

		fieldLength += covariateData.getDimension();

		addVariable(covariateData); // this can have missing values for imputation

	}

	//	@Override
	public void sSetupIntervals() {

		intervals.clear();
		intervals.ensureCapacity(fieldLength);

		
		NodeRef x;
		for (int i = 0; i < tree.getInternalNodeCount(); i++) {
			x = tree.getInternalNode(i);
			intervals.add(new CoalescentIntervalWithData(tree.getNodeHeight(x), Double.NaN, 0, CoalescentEventType.COALESCENT));
		}
		for (int i = 0; i < tree.getExternalNodeCount(); i++) {
			x = tree.getExternalNode(i);
			if (tree.getNodeHeight(x) > 1E-5){
				intervals.add(new CoalescentIntervalWithData(tree.getNodeHeight(x), Double.NaN, 0, CoalescentEventType.NEW_SAMPLE));
			}
		}
				
		dr.util.HeapSort.sort(intervals);
		for(int i = 0; i < intervals.size(); i++){
			intervals.get(i).lineage = getLineageCount(i);
		}
		
		
		
		for (int i = 0; i < covariateTimes.getDimension(); i++) {
			intervals.add(new CoalescentIntervalWithData(covariateTimes.getParameterValue(i),
					covariateData.getParameterValue(i), 0, CoalescentEventType.NOTHING));
		}
		dr.util.HeapSort.sort(intervals);
		
		
		
		double a = 0, b = 0;
		for (int i = 0; i < intervals.size(); i++) {
			b = intervals.get(i).length;
			intervals.get(i).length = intervals.get(i).length - a;
			a = b;
		}
		
				
		for(int i = 0; i < intervals.size(); i++){
			if(intervals.get(i).type.equals(CoalescentEventType.NOTHING)){
				int j = i - 1;
				double temp = intervals.get(i).datum;
				while(j > -1 && 
						!intervals.get(j).type.equals(CoalescentEventType.NOTHING)){
					intervals.get(j).datum = temp;
					j--;
				} 
			}
		}
		
		for(int i = 0; i < intervals.size(); i++){
			if(!intervals.get(i).type.equals(CoalescentEventType.NOTHING)){
				int lcount = intervals.get(i).lineage;
				int j = i - 1;
				while(j > -1 && intervals.get(j).type.equals(CoalescentEventType.NOTHING)){
					intervals.get(j).lineage = lcount;
					j--;
				}
			}
		}
		
		for(int i = 0; i < intervals.size(); i++){
			if(intervals.get(i).lineage == 0){
				intervals.get(i).lineage = 1;
			}
		}
			
		intervalsKnown = true;
	}
	
	public double calculateLogLikelihood(){
		
		double logLike = 0;
		
		for(CoalescentIntervalWithData interval : intervals){
			
			if(interval.lineage > 1){
								
				double lineageChoose2 = Binomial.choose2(interval.lineage);
						
				logLike += -lineageChoose2*Math.exp(-interval.datum)*interval.length;
				
				if(interval.type.equals(CoalescentEventType.COALESCENT)){
					logLike += -interval.datum;
				}
				
			}else{
				break;
			}
		}
		
		if(Double.isNaN(logLike)){
			System.out.println(logLike);
			System.out.println(intervals);
			System.out.println(tree.getNodeHeight(tree.getRoot()));
			System.exit(-1);
		}
		
		return logLike;
	}

	public void setupGMRFWeights() {
		super.setupGMRFWeights();
	}

	public void storeState() {
		storedIntervals = new ArrayList<CoalescentIntervalWithData>(intervals.size());
		for (CoalescentIntervalWithData interval : intervals) {
			storedIntervals.add(interval.clone());
		}
	}

	public void restoreState() {
		intervals = storedIntervals;
		storedIntervals.clear();
	}
	
	public String toString(){
		return intervals.toString();
	}
	
	public int getNumberOfIntervals(){
		return intervals.size();
	}
	
	public CoalescentIntervalWithData getDataInterval(int interval){
		return intervals.get(interval);
	}

	public class CoalescentIntervalWithData implements Comparable<CoalescentIntervalWithData>, Cloneable {
		public CoalescentEventType type;
		public double length;
		public int lineage;
		public double datum;

		public CoalescentIntervalWithData(double length, double datum, int lineage, CoalescentEventType type) {
			this.length = length;
			this.type = type;
			this.datum = datum;
			this.lineage = lineage;
		}

		public int compareTo(CoalescentIntervalWithData a) {

			if (a.length < this.length) {
				return 1;
			} else if (a.length == this.length) {
				Logger.getLogger("dr.evomodel.coalescent").severe("The current model " +
						"has 2 internal nodes or 1 node and 1 covariate at the same height\n" +
						a.toString() + "\n" + this.toString());
				return 0;
			}
			return -1;
		}

		public String toString() {
			return "(" + length + "," + type + "," + datum + "," + lineage + ")";
		}

		public CoalescentIntervalWithData clone() {
			return new CoalescentIntervalWithData(length, datum, lineage, type);
		}

	}
	
}
