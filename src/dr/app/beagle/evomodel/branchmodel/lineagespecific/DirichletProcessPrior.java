package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import java.util.ArrayList;
import java.util.List;


import dr.app.bss.Utils;
import dr.inference.model.Parameter;
import dr.math.distributions.Distribution;
import dr.math.distributions.MultivariateDistribution;
import dr.math.distributions.UniformDistribution;

public class DirichletProcessPrior implements MultivariateDistribution {

	private List<List<Parameter>> hyperparameterList;
	
	// distribution of realized parameter values
	private MultivariateDistribution distribution;
	// vector z of branch assignments (length 2n-1)
	private ArrayList<Integer> mapping;
	private Parameter categoriesParameter;
	// K
	private int categoryCount;
	// N=2n-1
	private int mappingLength;
	// concentration parameter of DPP
	private Parameter gamma;

	private Parameter uniquelyRealizedParameters;
	
	private double[] cachedLogFactors;
	
	public DirichletProcessPrior(MultivariateDistribution distribution, //
			ArrayList<Integer> mapping, //
			Parameter categoriesParameter,
			Parameter uniquelyRealizedParameters,
			List<List<Parameter>> hyperparameterList,
			Parameter gamma //
	) {

		
		this.hyperparameterList=hyperparameterList;
		
		this.uniquelyRealizedParameters = uniquelyRealizedParameters;
		
		this.distribution = distribution;
		this.mapping = mapping;

		this.categoryCount = Utils.findMaximum(mapping) + 1;// parameterValues.size();
		this.mappingLength = mapping.size();

		this.gamma = gamma;

		this.categoriesParameter = categoriesParameter;
		
		cacheLogFactorials();
		
	}// END: Constructor

	/**
	 * Assumes mappings start from index 0
	 * */
	private int[] getCounts() {

		// eta_k parameters (number of assignments to each category)
		int[] counts = new int[categoryCount];
		for (int i = 0; i < mappingLength; i++) {

			counts[mapping.get(i)]++;

		}// END: i loop

		return counts;
	}// END: getCounts

	public double getPriorLoglike() {

		int[] counts = getCounts();
		double loglike = categoryCount * Math.log(gamma.getParameterValue(0));

		for (int k = 0; k < categoryCount; k++) {

//			loglike += logfactor(counts[k] - 1);
			loglike+= cachedLogFactors[counts[k] - 1];
			

		}// END: k loop

		for (int i = 1; i <= mappingLength; i++) {
			loglike -= Math.log(gamma.getParameterValue(0) + i - 1);
		}// END: i loop

		return loglike;
	}// END: getPriorLoglike

	@Override
	public double logPdf(double[] x) {

		double loglike = 0.0;
		for (int i = 0; i < x.length; i++) {

//			loglike += distribution.logPdf(x[i]);

		}

		return loglike;
	}// END: logPdf

	
	private void cacheLogFactorials() {
		
        cachedLogFactors = new double[categoryCount];
        cachedLogFactors[0] = 0.0;
        for (int j = 1; j < categoryCount; j++) {
            cachedLogFactors[j] = cachedLogFactors[j - 1] + Math.log(j);
        }
		
	}//END: cacheFactorials
	
	private double logfactor(int n) {

//		double fact = 1;
//		   for ( int i = 1 ; i <= n ; i++ ) {
//	            fact = fact*i;
//	      }
		
		double logfactor = 0.0;

		for (int i = 1; i <= n; i++) {
			logfactor += Math.log(i);
		}

		return logfactor;
	}

	@Override
	public double[][] getScaleMatrix() {
		return null;
	}

	@Override
	public double[] getMean() {
		return null;
	}

	@Override
	public String getType() {
		return null;
	}

	public static ArrayList<Integer> chinRest(int N, double gamma) {

		ArrayList<Integer> mapping = new ArrayList<Integer>();

		mapping.add(0, 0);
		int nextFree = 1;

		for (int i = 1; i < N; i++) {

			double u1 = Math.random();
			if (u1 < (gamma / (gamma + i))) {

				// sit at new table
				mapping.add(i, nextFree);
				nextFree++;

			} else {

				// choose existing table with weights by number of customers
				int samplePos = -Integer.MAX_VALUE;
				double cumProb = 0.0;
				double u2 = Math.random();
				double size = (double) mapping.size();
				
				for(int j : mapping) {
					
					cumProb += 1/size;
					
					if (u2 < cumProb) {
						samplePos = j;
						break;
					}
					
				}
				
				mapping.add(i, samplePos);

			}// END: u1 check

		}// END: i loop

		return mapping;
	}// END: chinProc

//	public static void main(String args[]) {
//
//		// szybki lopez
//
//		int N = 10;
//		Parameter gamma = new Parameter.Default(0.9);
//		ArrayList<Integer> mapping = chinRest(N, gamma);
//		int K = Utils.findMaximum(mapping) + 1;
//
//		MultivariateDistribution dist =  null;// new UniformDistribution(0.0, 1.1);
//		
//		double[] parameterValues = new double[K];
//		for (int i = 0; i < K; i++) {
//
//			double u = Math.random();
//			parameterValues[i] = u;
//
//		}// END: i loop
//
//		DirichletProcessPrior dpp = new DirichletProcessPrior(dist,
//				mapping, gamma);
//
//		System.out.println("Mappings:");
//		Utils.printArray(mapping.toArray());
//		System.out.println("K=" + K);
//		// System.out.println("Values:");
//		// Utils.printArray(parameterValues);
//		System.out.println("Etas:");
//		Utils.printArray(dpp.getCounts());
//
//		double loglike = dpp.getPriorLoglike() +dpp.logPdf( parameterValues );
//
//		System.out.println("Loglikelihood=" + loglike);
//
//	}// END: main

}// END: class

