package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import java.util.ArrayList;

import dr.app.bss.Utils;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.math.GammaFunction;
import dr.math.distributions.Distribution;
import dr.math.distributions.LogNormalDistribution;
import dr.math.distributions.MultivariateDistribution;
import dr.math.distributions.UniformDistribution;

public class StickBreakingProcessPrior implements MultivariateDistribution {

	// distribution of realized parameter values
	private Distribution distribution;
	// vector z of branch assignments (length 2N-1)
	private ArrayList<Integer> mapping;
	// vector with realized parameter values (length K)
	private ArrayList<Double> parameterValues;

	// K
	private int categoryCount;
	// 2N-1
	private int mappingLength;
	// concentration parameter of DPP
	private double gamma;

	public StickBreakingProcessPrior(Distribution distribution, //
			ArrayList<Integer> mapping, //
			ArrayList<Double> parameterValues, //
			double gamma //
	) {

		this.distribution = distribution;
		this.mapping = mapping;
		this.parameterValues = parameterValues;

		this.categoryCount = parameterValues.size()+1;
		this.mappingLength = mapping.size();

		this.gamma = gamma;

	}// END: Constructor

	private int[] getCounts() {

		// eta_k parameters (number of assignments to each category)
		int[] counts = new int[categoryCount];
		for (int i = 0; i < mappingLength; i++) {

			counts[mapping.get(i)]++;

		}// END: i loop

		return counts;
	}// END: getCounts

	public double getPriorLoglike() {

		double loglike = categoryCount * Math.log(gamma);

		// eta_k parameters (number of assignments to each category)
		int[] counts = getCounts();
		for (int k = 0; k < categoryCount; k++) {

			loglike += Math.log(factorial(counts[k] - 1));

		}// END: k loop

		for (int i = 0; i < mappingLength; i++) {

			loglike -= Math.log(gamma + i - 1);

		}// END: i loop

		return loglike;
	}// END: getPriorLoglike

	@Override
	public double logPdf(double[] x) {

		double loglike = 0.0;

		for (int i = 0; i < categoryCount; i++) {

			
			loglike += distribution.logPdf(parameterValues.get(i));
			
//			System.err.println("FUBAR" + loglike);
			
		}

		return loglike;
	}// END: logPdf

	private double factorial(double n) {

		if (n <= 0.0) {
			return 1.0;
		} else {
			return n * factorial(n - 1);
		}

	}// END: factorial

	public double getMultivariateLoglike() {

		double logLike = 0.0;

		int[] counts = getCounts();
		double countSum = Utils.sumArray(counts);

		logLike += GammaFunction.lnGamma(countSum);
		for (int i = 0; i < categoryCount; i++) {

			logLike += ((counts[i] - 1) * Math.log(parameterValues.get(i)) - GammaFunction
					.lnGamma(counts[i]));

		}// END: i loop

		return logLike;
	}// END: getMultivariateLoglike
	
	// S(n,k) = S(n-1,k-1) * k*S(n-1, k)
	public int stirlingNumberSecondKind(int n, int k) {

		if (k == 1 || k == n) {
			return 1;
		}

		return stirlingNumberSecondKind(n - 1, k - 1) + k
				* stirlingNumberSecondKind(n - 1, k);
	}// END: stirlingNumberSecondKind

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

		ArrayList<Integer> mapping = new ArrayList<>();

		mapping.add(0, 1);
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

				for (int j = 0; j < mapping.size(); j++) {

					cumProb += (1 / (j + 1));

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

	public static void main(String args[]) {

		// szybki lopez

		int N = 10;
		double gamma = 2.9;
		ArrayList<Integer> mappings = chinRest(N, gamma);
		int K = Utils.findMaximum(mappings);

		Distribution dist = new UniformDistribution(0.0, 1.1);

		ArrayList<Double> parameterValues = new ArrayList<>();
		for (int i = 0; i < K; i++) {

			double u = Math.random();
			parameterValues.add(i, u);

		}// END: i loop

		StickBreakingProcessPrior dpp = new StickBreakingProcessPrior(dist, mappings, parameterValues, gamma);

		
		System.out.println("Mappings:");
		Utils.printArray(mappings.toArray());
		System.out.println("K=" + K);
		System.out.println("Values:");
		Utils.printArray(parameterValues.toArray());
		System.out.println("Etas:");
		Utils.printArray(dpp.getCounts());
		
		double loglike =  dpp.getPriorLoglike(); //+dpp.logPdf( null );

		System.out.println("Loglikelihood=" + loglike);	
		
	}// END: main

}// END: class

