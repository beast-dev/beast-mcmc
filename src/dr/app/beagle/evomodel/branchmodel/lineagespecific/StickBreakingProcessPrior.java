package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import java.util.ArrayList;

import dr.app.bss.Utils;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.math.GammaFunction;
import dr.math.distributions.Distribution;
import dr.math.distributions.MultivariateDistribution;

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

		this.categoryCount = parameterValues.size();
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

	public int stirlingNumberFirstKind(int n, int k) {

		if (k == 1 || k == n) {
			return 1;
		}

		return stirlingNumberFirstKind(n - 1, k - 1) + k
				* stirlingNumberFirstKind(n - 1, k);
	}// END: stirlingNumberFirstKind

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

}// END: class

