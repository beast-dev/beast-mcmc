package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import org.apache.commons.math.MathException;

import dr.app.bss.Utils;
import dr.math.MathUtils;
import dr.math.UnivariateFunction;
import dr.math.distributions.BetaDistribution;
import dr.math.distributions.Distribution;

public class StickBreakingProcessPrior implements Distribution {

	private double intensity;
	private BetaDistribution beta;
	private int K;

	public StickBreakingProcessPrior(double intensity, int K) {

		this.intensity = intensity;
		this.K = K;

		beta = new BetaDistribution(this.intensity, 1.0);

	}// END: Constructor

	public double[] stickBreak() throws MathException {

		double[] probs = new double[K];
		double stickLength = 1.0;
		for (int i = 0; i < K; i++) {

			double beta = rBeta();
			probs[i] = beta * stickLength;
			stickLength *= (1 - beta);

		}

		return probs;
	}

	@Override
	public double pdf(double x) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double logPdf(double x) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double cdf(double x) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double quantile(double y) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double mean() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double variance() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public UnivariateFunction getProbabilityDensityFunction() {
		// TODO Auto-generated method stub
		return null;
	}

	private double rBeta() throws MathException {
		return beta.inverseCumulativeProbability(MathUtils.nextDouble());
	}// END: rbeta

	public static void main(String args[]) {

		try {

			int K = 10;
			double intensity = 1.0;

			StickBreakingProcessPrior sbp = new StickBreakingProcessPrior(
					intensity, K);
			double[] probs = sbp.stickBreak();

			System.out.println(Utils.sumArray(probs));

		} catch (MathException e) {
			e.printStackTrace();
		}

	}// END: main

}// END: class
