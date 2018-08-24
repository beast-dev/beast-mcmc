/*
 * DirichletProcessOperator.java
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

package dr.evomodel.branchmodel.lineagespecific;

import org.apache.commons.math.MathException;

import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

@SuppressWarnings("serial")
public class DirichletProcessOperator extends SimpleMCMCOperator implements
		GibbsOperator {

	private static final boolean DEBUG = false;

	private DirichletProcessPrior dpp;

	private int realizationCount;
	private int uniqueRealizationCount;
	private double intensity;
	private int mhSteps;

	private Parameter categoriesParameter;
	 private CountableRealizationsParameter allParameters;
	private Parameter uniqueParameters;

	private CompoundLikelihood likelihood;

	public DirichletProcessOperator(DirichletProcessPrior dpp, //
			Parameter categoriesParameter, //
			Parameter uniqueParameters, //
			 CountableRealizationsParameter allParameters,
			Likelihood likelihood, //
			int mhSteps, //
			double weight//
	) {

		this.dpp = dpp;
		this.intensity = dpp.getGamma();
		this.uniqueRealizationCount = dpp.getCategoryCount();
		this.realizationCount = categoriesParameter.getDimension();

		this.categoriesParameter = categoriesParameter;
		 this.allParameters = allParameters;
		this.uniqueParameters = uniqueParameters;
		this.likelihood = (CompoundLikelihood) likelihood;
//		this.likelihood =  likelihood;
		
		this.mhSteps = mhSteps;
		
		setWeight(weight);

	}// END: Constructor

	public Parameter getParameter() {
		return categoriesParameter;
	}// END: getParameter

	public Variable getVariable() {
		return categoriesParameter;
	}// END: getVariable

	@Override
	public double doOperation() {

		try {

//			doOperate();
			doOp();

		} catch (MathException e) {
			e.printStackTrace();
		}// END: try-catch block

		return 0.0;
	}// END: doOperation

	private void doOp() throws MathException {
		
		for (int index = 0; index < realizationCount; index++) {
		
			int[] occupancy = new int[uniqueRealizationCount];
			for (int i = 0; i < realizationCount; i++) {
				if (i != index) {
					int j = (int) categoriesParameter.getParameterValue(i);
					occupancy[j]++;
				}// END: i check
			}// END: i loop

			
	        double[] existingValues = new double[uniqueRealizationCount];
	        int counter = 0;
	        int singletonIndex = -1;
	        for(int i = 0; i < uniqueRealizationCount;i++){
	            if(occupancy[i] > 0) {
	            	
	                occupancy[counter] = occupancy[i];
	                existingValues[counter++] = dpp.getUniqueParameter(i) .getParameterValue(0);

	            } else {
	            
	            	singletonIndex = i;

	            }//END: occupancy check
	            
	        }//END: i loop
			
			
			// Propose new value(s)
			double[] baseProposals = new double[realizationCount];
			for (int i = 0; i < baseProposals.length; i++) {
				
				baseProposals[i] = dpp.baseModel.nextRandom()[0];
				
			}
			
	        // If a singleton
            if(singletonIndex > -1) {
            	
                baseProposals[0] = uniqueParameters.getParameterValue(singletonIndex);

            }

			double[] logClusterProbs = new double[uniqueRealizationCount];
            
            // draw existing
            int i;
            for(i = 0; i < counter; i++) {
            
            	  logClusterProbs[i] = Math.log(occupancy[i] / (realizationCount - 1 + intensity));
            	  
            	  double value =  allParameters.getParameterValue(index);
            	  double candidate = existingValues[i];
            	  allParameters.setParameterValue(index, candidate);
				  likelihood.makeDirty();
            	  
            	  logClusterProbs[i] = logClusterProbs[i] + likelihood.getLikelihood(index) .getLogLikelihood();
//				  logClusterProbs[i] = logClusterProbs[i] + likelihood .getLogLikelihood();
				  
//            	  System.out.println(likelihood.getLikelihood(index) .getLogLikelihood() + " " + likelihood .getLogLikelihood());
            	  
            	  allParameters.setParameterValue(index, value);
            	  likelihood.makeDirty();
            	  
            }
            
            // draw new
            for(; i < logClusterProbs.length; i++){

            	logClusterProbs[i] = Math.log((intensity) / (realizationCount - 1 + intensity)); 
//            	logClusterProbs[i] = Math.log(intensity / uniqueRealizationCount / (realizationCount - 1 + intensity));
            	
            	  double value =  allParameters.getParameterValue(index);
            	 double candidate = baseProposals[i - counter];
            	 allParameters.setParameterValue(index, candidate);
 				 likelihood.makeDirty();
 				 
            	 logClusterProbs[i] = logClusterProbs[i] + likelihood.getLikelihood(index).getLogLikelihood();
// 				 logClusterProbs[i] = logClusterProbs[i] + likelihood.getLogLikelihood();
 				 
//           	  System.out.println(likelihood.getLikelihood(index) .getLogLikelihood() + " " + likelihood .getLogLikelihood());
            	 
            	 allParameters.setParameterValue(index, value);
            	 likelihood.makeDirty();
            	 
            }
            
            double smallestVal = logClusterProbs[0];
            for(i = 1; i < uniqueRealizationCount; i++){
                
            	if(smallestVal > logClusterProbs[i]) {
                    smallestVal = logClusterProbs[i];
                }
            
            }
            
            
            double[] clusterProbs = new double[uniqueRealizationCount];
            for(i = 0; i < clusterProbs.length;i++) {
                    clusterProbs[i] = Math.exp(logClusterProbs[i]-smallestVal);
            
            }

//            dr.app.bss.Utils.printArray(clusterProbs);
//         	System.exit(-1);
            
			// sample
			int sampledCluster = MathUtils.randomChoicePDF(clusterProbs);
			categoriesParameter.setParameterValue(index, sampledCluster);
            
            
		}//END: index loop
		
		
	}//END: doOp
	
	
	private void doOperate() throws MathException {

		// int index = 0;
		for (int index = 0; index < realizationCount; index++) {

			int[] occupancy = new int[uniqueRealizationCount];
			for (int i = 0; i < realizationCount; i++) {

				if (i != index) {

					int j = (int) categoriesParameter.getParameterValue(i);
					occupancy[j]++;

				}// END: i check

			}// END: i loop

			if (DEBUG) {
				System.out.println("N[-index]: ");
				dr.app.bss.Utils.printArray(occupancy);
			}

			Likelihood clusterLikelihood = (Likelihood) likelihood.getLikelihood(index);
//			Likelihood clusterLikelihood = likelihood;
			
			int category = (int) categoriesParameter.getParameterValue(index);
			double value = uniqueParameters.getParameterValue(category);
			
			double[] clusterProbs = new double[uniqueRealizationCount];
			
			for (int i = 0; i < uniqueRealizationCount; i++) {

				double logprob = 0;
				if (occupancy[i] == 0) {// draw new

					// draw from base model, evaluate at likelihood

					double candidate = dpp.baseModel.nextRandom()[0];

					uniqueParameters.setParameterValue(category, candidate);
					double loglike = clusterLikelihood.getLogLikelihood();
					uniqueParameters.setParameterValue(category, value);

					logprob = Math.log((intensity) / (realizationCount - 1 + intensity)) + loglike;

				} else {// draw existing

					// likelihood for component x_index

					double candidate = dpp.getUniqueParameter(i)
							.getParameterValue(0);

					uniqueParameters.setParameterValue(category, candidate);
					double loglike = clusterLikelihood.getLogLikelihood();
					uniqueParameters.setParameterValue(category, value);

					logprob = Math.log(occupancy[i]) / (realizationCount - 1 + intensity) + loglike;

				}// END: occupancy check

				clusterProbs[i] = logprob;
			}// END: i loop

			dr.app.bss.Utils.exponentiate(clusterProbs);

			if (DEBUG) {
				System.out.println("P(z[index] | z[-index]): ");
				dr.app.bss.Utils.printArray(clusterProbs);
			}

			// sample
			int sampledCluster = MathUtils.randomChoicePDF(clusterProbs);
			categoriesParameter.setParameterValue(index, sampledCluster);

			if (DEBUG) {
				System.out
						.println("sampled category: " + sampledCluster + "\n");
			}

		}// END: index loop

	}// END: doOperate

	@Override
	public String getOperatorName() {
		return DirichletProcessOperatorParser.DIRICHLET_PROCESS_OPERATOR;
	}

	@Override
	public String getPerformanceSuggestion() {
		return null;
	}// END: getPerformanceSuggestion
	
}// END: class
