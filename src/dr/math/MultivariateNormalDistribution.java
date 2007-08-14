package dr.math;

import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.SymmetricMatrix;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Jun 13, 2007
 * Time: 1:47:26 PM
 * To change this template use File | Settings | File Templates.
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
        logDet = calculatePrecisionMatrixDeterminate(precision);

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
            return (new Matrix(precision).determinant());
        } catch (
                IllegalDimension e) {
            throw new RuntimeException(e.getMessage());
        }
    }

    public double logPdf(Parameter x) {
        return logPdf(x.getParameterValues(), mean, precision, logDet, 1.0);
    }

    public static final double logPdf(double[] x, double[] mean, double[][] precision,
                                      double logDet, double scale) {
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

//	public double[][] inverseScaleMatrix() {
//			return inverseScaleMatrix;
//		}

    public double[] nextMultivariateNormal() {
        return nextMultivariateNormal(mean, precision);
    }


    public static double[] nextMultivariateNormal(double[] mean, double[][] precision) {

        int dim = mean.length;

        double[][] cholesky = null;

        double[][] variance = new SymmetricMatrix(precision).inverse().toComponents();

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
