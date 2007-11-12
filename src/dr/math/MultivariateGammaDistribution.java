package dr.math;

import dr.inference.model.Parameter;

/**
 * @author Marc Suchard
 */
public class MultivariateGammaDistribution implements MultivariateDistribution {

	// todo Currently this implements a product of independent Gammas, need to re-code as true multivariate distribution

	public static final String TYPE = "multivariateGamma";

	public MultivariateGammaDistribution(double[] shape, double[] scale) {

		if (shape.length != scale.length)
			throw new RuntimeException("Creation error in MultivariateGammaDistribution");

		dim = shape.length;

		this.shape = shape;
		this.scale = scale;

	}

	public double logPdf(Parameter x) {

		double logPdf = 0;

		for (int i = 0; i < dim; i++)
			logPdf += GammaDistribution.logPdf(x.getParameterValue(i),
					shape[i], scale[i]);

		return logPdf;
	}

	public double[][] getScaleMatrix() {
		throw new RuntimeException("Not yet implemented");
	}

	public double[] getMean() {
		throw new RuntimeException("Not yet implemented");
	}

	public String getType() {
		return TYPE;
	}

	private double[] shape;
	private double[] scale;
	private int dim;
}
