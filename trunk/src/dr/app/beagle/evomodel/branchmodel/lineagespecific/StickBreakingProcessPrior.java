package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import org.apache.commons.math.MathException;

import dr.app.bss.Utils;
import dr.math.MathUtils;
import dr.math.UnivariateFunction;
import dr.math.distributions.BetaDistribution;
import dr.math.distributions.Distribution;
import dr.util.HeapSort;

public class StickBreakingProcessPrior implements Distribution {

	private double intensity;
	private BetaDistribution beta;
	private int K;
	private double[] probabilities;
	private int[] support;
	private UnivariateFunction pdfFunction;

	public StickBreakingProcessPrior(double intensity, int K) {

		this.intensity = intensity;
		this.K = K;

		beta = new BetaDistribution(this.intensity, 1.0);

	}// END: Constructor

	public void stickBreak() throws MathException {

		probabilities = new double[K];
		support = new int[K];

		double stickLength = 1.0;
		int i = 0;
		for (int x = 1; x <= K; x++) {

			support[i] = x;

			double beta = rBeta();
			probabilities[i] = beta * stickLength;
			stickLength *= (1 - beta);

			i = i + 1;

		}

		pdfFunction = new UnivariateFunction() {
			public double evaluate(double x) {
				return pdf(x);
			}

			public double getLowerBound() {
				return 1.0;
			}

			public double getUpperBound() {
				return K;
			}
		};

	}// END: stickBreak

	private boolean inSupport(double x) {

		for (int i = 0; i < K; i++) {

			if (support[i] == x) {
				return true;
			}

		}

		return false;
	}// END: inSupport

	// probability mass function
	@Override
	public double pdf(double x) {

		if (inSupport(x)) {

			return probabilities[(int) (x - 1)];

		}

		return 0;
	}// END: pdf

	@Override
	public double logPdf(double x) {
		return Math.log(pdf(x));
	}

	@Override
	public double cdf(double x) {

		double sum = 0.0;
		for (int i = 0; i < K; i++) {

			if (i <= x) {
				sum += probabilities[i];
			}

		}

		return sum;
	}// END: cdf

	@Override
	public double quantile(double y) {

		if (y < 0.0 || y > 1.0) {
			throw new RuntimeException("Quantile out of range.");
		}

		// if (y == 0.0) {
		// return support[indices[0]] - 1.0;
		// }

		int[] indices = new int[K];
		HeapSort.sort(probabilities, indices);

		// Utils.printArray(indices);

		return support[indices[(int) Math.ceil(y * K) - 1]];
	}

	@Override
	public double mean() {

		double mean = 0.0;
		for (int i = 0; i < K; i++) {

			mean += (support[i] * probabilities[i]);

		}

		return mean;
	}// END: mean

	@Override
	public double variance() {

		double mean = mean();

		double variance = 0.0;
		for (int i = 0; i < K; i++) {

			variance += Math.pow((support[i] - mean), 2) * probabilities[i];

		}

		return variance;
	}// END: variance

	@Override
	public UnivariateFunction getProbabilityDensityFunction() {
		return pdfFunction;
	}

	private double rBeta() throws MathException {
		return beta.inverseCumulativeProbability(MathUtils.nextDouble());
	}// END: rbeta

	public double[] getProbabilities() {
		return probabilities;
	}

	public int[] getSupport() {
		return support;
	}

	public static void main(String args[]) {

		try {

			int K = 10; // close to infinity but not too big to save compute
						// time
			double intensity = 2;

			StickBreakingProcessPrior sbp = new StickBreakingProcessPrior(
					intensity, K);
			sbp.stickBreak();

			double[] probs = sbp.getProbabilities();

			Utils.print2Arrays(sbp.getSupport(), probs, probs.length);

			// System.out.println();
			// System.out.println(Utils.sumArray(probs));

			System.out.println();

			System.out.println(sbp.quantile(0.5));

		} catch (Exception e) {
			e.printStackTrace();
		}

	}// END: main

}// END: class
