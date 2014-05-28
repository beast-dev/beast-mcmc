package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import dr.app.bss.Utils;
import dr.inference.model.Parameter;
import dr.math.GammaFunction;
import dr.math.distributions.Distribution;
import dr.math.distributions.MultivariateDistribution;

public class StickBreakingProcessPrior implements MultivariateDistribution {

	private Distribution baseDistribution;
	private Parameter mdParameter;
	private Parameter uCategories;
	private int categoryCount;

	public StickBreakingProcessPrior(Distribution baseDistribution,//
			Parameter mdParameter, //
			Parameter uCategories //
	) {

		this.baseDistribution = baseDistribution;
		this.mdParameter = mdParameter;
		this.uCategories = uCategories;
		this.categoryCount = mdParameter.getDimension();

	}// END: Constructor

	private int[] getCounts() {

		int[] counts = new int[categoryCount];
		int[] branchAssignments = new int[uCategories.getDimension()];

		for (int i = 0; i < branchAssignments.length; i++) {
			counts[branchAssignments[i]]++;
		}

		return counts;
	}// END: getBranchAssignmentCounts

	private double[] getValues() {

		double[] values = new double[categoryCount];

		for (int i = 0; i < categoryCount; i++) {

			values[i] = mdParameter.getParameterValue(i);

		}

		return values;
	}// END: getValues

	/*
	 * Distribution Likelihood
	 */
	public double getLogLikelihood() {

		double logLike = 0.0;

		int[] counts = getCounts();
		double[] values = getValues();
		double countSum = Utils.sumArray(counts);

		logLike += GammaFunction.lnGamma(countSum);
		for (int i = 0; i < categoryCount; i++) {

			logLike += ((counts[i] - 1) * Math.log(values[i]) - GammaFunction
					.lnGamma(counts[i]));

		}// END: i loop

		return logLike;
	}// END: getLogLikelihood

	@Override
	public double logPdf(double[] x) {
		double logLike = 0.0;

		int[] counts = getCounts();
		double countSum = Utils.sumArray(counts);

		logLike += GammaFunction.lnGamma(countSum);
		for (int i = 0; i < categoryCount; i++) {

			logLike += ((counts[i] - 1) * Math.log(x[i]) - GammaFunction
					.lnGamma(counts[i]));

		}// END: i loop

		return logLike;
	}// END: logPdf

	@Override
	public double[] getMean() {

		int[] counts = getCounts();
		double countSum = Utils.sumArray(counts);

		double[] mean = new double[categoryCount];

		for (int i = 0; i < categoryCount; i++) {
			mean[i] = counts[i] / countSum;
		}

		return mean;
	}// END: mean

	@Override
	public double[][] getScaleMatrix() {
		return null;
	}

	@Override
	public String getType() {
		return null;
	}

}// END: class

