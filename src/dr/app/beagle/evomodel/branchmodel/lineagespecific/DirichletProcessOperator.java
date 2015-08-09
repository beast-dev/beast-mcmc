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

package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import org.apache.commons.math.MathException;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;


/**
 * @author Filip Bielejec
 */
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
//	private CountableRealizationsParameter countableRealizationsParameter;
	
	// Compound Parameter 
	private CompoundParameter parameter;
	private CompoundLikelihood likelihood;
	
	public DirichletProcessOperator(DirichletProcessPrior dpp, 
			Parameter categoriesParameter, 
//			CountableRealizationsParameter countableRealizationsParameter,
			CompoundParameter parameter,
		    CompoundLikelihood likelihood,
			int mhSteps,
			double weight) {

		this.dpp = dpp;
		this.intensity = dpp.getGamma();
		this.uniqueRealizationCount = dpp.getCategoryCount();
		this.realizationCount = categoriesParameter.getDimension();
		
		this.categoriesParameter = categoriesParameter;
//		this.countableRealizationsParameter = countableRealizationsParameter;
		this.parameter = parameter;
		this.likelihood = likelihood;
		
		this.mhSteps = mhSteps;
		setWeight(weight);

	}// END: Constructor

    public Parameter getParameter() {
        return categoriesParameter;
    }//END: getParameter
	
    public Variable<?> getVariable() {
        return categoriesParameter;
    }//END: getVariable
	
	@Override
	public double doOperation() throws OperatorFailedException {

		try {

			doOperate();

		} catch (MathException e) {
			e.printStackTrace();
		}// END: try-catch block

		return 0.0;
	}// END: doOperation

	private void doOperate() throws MathException {

//		System.out.println(likelihood.toString());
//		System.out.println("likelihood count: " + likelihood.getLikelihoodCount());
		
//		int index = 0;
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
			
			double[] clusterProbs = new double[uniqueRealizationCount];
			for (int i = 0; i < uniqueRealizationCount; i++) {

				double loglike = 0;
				double logprob = 0;
				if (occupancy[i] == 0) {// draw new

					// M-H for poor people
					// draw from base model, evaluate at likelihood
					double candidate = 0.0;
					for (int j = 0; j < mhSteps; j++) {
						 candidate = dpp.baseModel.nextRandom()[0];
						 loglike += Math.exp(getPartialLoglike(index, candidate));
					}
					
					loglike /= mhSteps;
					logprob = Math.log( ( (intensity) / (realizationCount - 1 + intensity) ) * loglike   ) ;

				} else {// draw existing
					
					// likelihood for component x_index
					 double value = dpp.getUniqueParameter(i).getParameterValue(0);
					 loglike = getPartialLoglike(index, value);
					
					 logprob = Math.log(occupancy[i]) / (realizationCount - 1 + intensity) + loglike;

				}// END: occupancy check
				
				clusterProbs[i] = logprob;
			}// END: i loop

			//rescale
			double max = dr.app.bss.Utils.max(clusterProbs);
			for (int i = 0; i < clusterProbs.length; i++) {
				clusterProbs[i] -=  max;
			}

			dr.app.bss.Utils.exponentiate(clusterProbs);
//			dr.app.bss.Utils.normalize(clusterProbs);
			
			if (DEBUG) {
				System.out.println("P(z[index] | z[-index]): ");
				dr.app.bss.Utils.printArray(clusterProbs);
			}

			// sample
			int sampledCluster = MathUtils.randomChoicePDF(clusterProbs);
			categoriesParameter.setParameterValue(index, sampledCluster);

			if (DEBUG) {
				System.out.println("sampled category: " + sampledCluster + "\n");
			}

		}// END: index loop

	}// END: doOperate

	private double getPartialLoglike(int index, double candidate) {

		DistributionLikelihood dl = (DistributionLikelihood) likelihood .getLikelihood(index);
		
		int category = (int) categoriesParameter.getParameterValue(index);
		double value = parameter.getParameterValue(category);
		
		double loglike = 0.0;
		if (candidate != value) {

			//TODO: which category
			parameter.setParameterValue(category, candidate);
			loglike = dl.getLogLikelihood();
			parameter.setParameterValue(category, value);

		} else {

			loglike = dl.getLogLikelihood();

		}

		return loglike;
	}// END: getPartialLoglike
	
	@Override
	public String getOperatorName() {
		return DirichletProcessOperatorParser.DIRICHLET_PROCESS_OPERATOR;
	}

	@Override
	public String getPerformanceSuggestion() {
		return null;
	}// END: getPerformanceSuggestion

	@Override
	public int getStepCount() {
		return realizationCount;
	}// END: getStepCount

}// END: class
