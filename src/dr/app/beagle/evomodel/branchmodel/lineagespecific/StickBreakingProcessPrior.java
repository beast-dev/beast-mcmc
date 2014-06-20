package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import java.util.ArrayList;

import dr.app.bss.Utils;
import dr.math.distributions.Distribution;
import dr.math.distributions.MultivariateDistribution;
import dr.math.distributions.UniformDistribution;

public class StickBreakingProcessPrior implements MultivariateDistribution {

	// distribution of realized parameter values
	private Distribution distribution;
	// vector z of branch assignments (length 2N-1)
	private ArrayList<Integer> mapping;

	// K
	private int categoryCount;
	// N=2n-1
	private int mappingLength;
	// concentration parameter of DPP
	private double gamma;

	public StickBreakingProcessPrior(Distribution distribution, //
			ArrayList<Integer> mapping, //
			double gamma //
	) {

		this.distribution = distribution;
		this.mapping = mapping;

		this.categoryCount = Utils.findMaximum(mapping) + 1;// parameterValues.size();
		this.mappingLength = mapping.size();

		this.gamma = gamma;

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

		// eta_k parameters (number of assignments to each category)
		int[] counts = getCounts();

		double loglike = categoryCount * Math.log(gamma);

		for (int k = 0; k < categoryCount; k++) {

			loglike += logfactorial(counts[k] - 1);

		}// END: k loop

		for (int i = 1; i <= mappingLength; i++) {
			loglike -= Math.log(gamma + i - 1);
		}// END: i loop

		return loglike;
	}// END: getPriorLoglike

	@Override
	public double logPdf(double[] x) {

		double loglike = 0.0;
		for (int i = 0; i < x.length; i++) {

			loglike += distribution.logPdf(x[i]);

		}

		return loglike;
	}// END: logPdf

	private double logfactorial(int n) {

		double logfactorial = 0.0;

		for (int j = 1; j < n; j++) {
			logfactorial += +Math.log(j);
		}

		return logfactorial;
	}

	// TODO: test and remove
	private long factorial(long n) {

		if (n <= 0) {
			return 1;
		} else {
			return n * factorial(n - 1);
		}

	}// END: factorial

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

		int N = 3;
		double gamma = 2.9;
		ArrayList<Integer> mapping = chinRest(N, gamma);
		int K = Utils.findMaximum(mapping) + 1;

		Distribution dist = new UniformDistribution(0.0, 1.1);

		double[] parameterValues = new double[K];
		for (int i = 0; i < K; i++) {

			double u = Math.random();
			parameterValues[i] = u;

		}// END: i loop

		StickBreakingProcessPrior dpp = new StickBreakingProcessPrior(dist,
				mapping, gamma);

		System.out.println("Mappings:");
		Utils.printArray(mapping.toArray());
		System.out.println("K=" + K);
		// System.out.println("Values:");
		// Utils.printArray(parameterValues);
		System.out.println("Etas:");
		Utils.printArray(dpp.getCounts());

		double loglike = dpp.getPriorLoglike();// +dpp.logPdf( null );

		System.out.println("Loglikelihood=" + loglike);

	}// END: main

}// END: class

