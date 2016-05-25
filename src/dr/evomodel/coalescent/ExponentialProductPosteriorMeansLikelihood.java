/*
 * ExponentialProductPosteriorMeansLikelihood.java
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

import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;

/**
 * Calculates a product of exponential densities and exponential tail probabilities.
 *
 * @author Guy Baele
 */

public class ExponentialProductPosteriorMeansLikelihood extends Likelihood.Abstract {
	
	//not used at the moment
	public final static boolean FIXED_TREE = false;
    public static final boolean DEBUG = false;
	
	private TreeModel treeModel;
	private double[] posteriorMeans;

	//make sure to use in combination with coalescentEventsStatistic
	public ExponentialProductPosteriorMeansLikelihood(TreeModel treeModel, double[] posteriorMeans) {
		super(treeModel);
		this.treeModel = treeModel;
		this.posteriorMeans = posteriorMeans;
	}
	
	public double calculateLogLikelihood() {
		
		//System.err.println(treeModel);
		
		double logPDF = 0.0;
		
		//System.err.println("log(popSize) = " + this.popSize);
		
		CoalescentTreeIntervalStatistic ctis = new CoalescentTreeIntervalStatistic(treeModel);
		int coalescentIntervalCounter = 0;

        if (DEBUG) {
            System.err.println("ExponentialProductPosteriorMeansLikelihood dimension: " + ctis.getDimension());
        }
		
		for (int i = 0; i < ctis.getDimension(); i++) {
			
			int combinations = (int)ctis.getLineageCount(i)*((int)ctis.getLineageCount(i)-1)/2;
			double branchLength = ctis.getStatisticValue(i);
			
			//System.err.println("combinations = " + combinations);
			//System.err.println("branchLength = " + branchLength);
			
			//System.err.println(ctis.getLineageCount(i));
			
			//single-lineage intervals are not counted
			if (ctis.getLineageCount(i) != 1) {
				
				//System.err.println(i + " -> lineage count: " + ctis.getLineageCount(i));
				
				if (i == ctis.getDimension()-1) {
					//coalescent event at root: exponential density
					//System.err.print("coalescent event at root: ");
					double logContribution = -posteriorMeans[coalescentIntervalCounter] - combinations*branchLength*Math.exp(-posteriorMeans[coalescentIntervalCounter]);
					if (DEBUG) {
                        System.err.println(i + ": " + logContribution);
                    }
                    logPDF += logContribution;
					//System.err.println(logContribution);
				} else if (ctis.getLineageCount(i) > ctis.getLineageCount(i+1)) {
					//coalescent event: exponential density
					//System.err.print("coalescent event (not at root): ");
					double logContribution = -posteriorMeans[coalescentIntervalCounter] - combinations*branchLength*Math.exp(-posteriorMeans[coalescentIntervalCounter]);
                    if (DEBUG) {
                        System.err.println(i + ": " + logContribution);
                    }
                    logPDF += logContribution;
					//System.err.println(logContribution);
					coalescentIntervalCounter++;
				} else {
					//sampling event: exponential tail probability
					//System.err.print("sampling event: ");
					double logContribution = -combinations*branchLength*Math.exp(-posteriorMeans[coalescentIntervalCounter]);
                    if (DEBUG) {
                        System.err.println(i + ": " + logContribution);
                    }
                    logPDF += logContribution;
					//System.err.println(logContribution);
				}
				
			}
			
		}
		
		//System.err.println("expoLike = " + logPDF + "\n");
		
		return logPDF;
		
	}

	/**
	 * Overridden to always return false.
	 */
	protected boolean getLikelihoodKnown() {
		return false;
	}

}
