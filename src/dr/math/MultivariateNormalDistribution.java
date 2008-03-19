package dr.math;

import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.SymmetricMatrix;

/**
 * @author Marc Suchard
 */
public class MultivariateNormalDistribution implements MultivariateDistribution {

	public static final String TYPE = "MultivariateNormal";

	private double[] mean;
	private double[][] precision;
	//	private int dim;
	private double logDet;

	public MultivariateNormalDistribution(double[] mean, double[][] precision) {
		this.mean = mean;
		this.precision = precision;
		logDet = Math.log(calculatePrecisionMatrixDeterminate(precision));
	}

	public String getType() {
		return TYPE;
	}

	public double[][] getScaleMatrix() {
		return precision;
	}

	public double[] getMean() {
		return mean;
	}

	public static final double calculatePrecisionMatrixDeterminate(double[][] precision) {
		try {
//            double det = new Matrix(precision).determinant();
//            if( det < 0 ) {
//                System.err.println("Negative determinant.  how?");
//                System.exit(-1);
//            }
//            return //Math.log(new Matrix(precision).determinant());
//                det;
			return new Matrix(precision).determinant();
		} catch (
				IllegalDimension e) {
			throw new RuntimeException(e.getMessage());
		}
	}

	public double logPdf(double[] x) {
		return logPdf(x, mean, precision, logDet, 1.0);
	}

	public double logPdf(Parameter x) {
		return logPdf(x.getParameterValues(), mean, precision, logDet, 1.0);
	}

	public static final double logPdf(double[] x, double[] mean, double[][] precision,
	                                  double logDet, double scale) {

		if (logDet == Double.NEGATIVE_INFINITY)
			return logDet;

		int dim = x.length;
		double[] delta = new double[dim];
		double[] tmp = new double[dim];

		for (int i = 0; i < dim; i++) {
			delta[i] = x[i] - mean[i];
		}

		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				tmp[i] += delta[j] * precision[j][i];
			}
		}

		double SSE = 0;

		for (int i = 0; i < dim; i++)
			SSE += tmp[i] * delta[i];

		SSE /= scale;

		return dim * logNormalize + 0.5 * logDet - 0.5 * SSE;

	}

	public double[] nextMultivariateNormal() {
		return nextMultivariateNormalPrecision(mean, precision);
	}


	public static double[] nextMultivariateNormalPrecision(double[] mean, double[][] precision) {

		double[][] variance = new SymmetricMatrix(precision).inverse().toComponents();
		return nextMultivariateNormalVariance(mean, variance);

	}

	public static double[] nextMultivariateNormalVariance(double[] mean, double[][] variance) {

		int dim = mean.length;

		double[][] cholesky = null;

		try {
			cholesky = (new CholeskyDecomposition(variance)).getL();
		} catch (IllegalDimension illegalDimension) {
			// todo - check for square variance matrix before here
		}

		double[] x = new double[dim];
		for (int i = 0; i < dim; i++)
			x[i] = mean[i];

		double[] epsilon = new double[dim];
		for (int i = 0; i < dim; i++)
			epsilon[i] = MathUtils.nextGaussian();

		for (int i = 0; i < dim; i++) {
			for (int j = i; j < dim; j++) {
				x[i] += cholesky[j][i] * epsilon[j];
				// caution: decomposition returns lower triangular
			}
		}
		return x;
	}

	public static final double logNormalize = -0.5 * Math.log(2.0 * Math.PI);

}
