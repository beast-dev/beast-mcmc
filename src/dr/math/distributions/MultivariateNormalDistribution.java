package dr.math.distributions;

import dr.math.MathUtils;
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
            return new Matrix(precision).determinant();
        } catch (IllegalDimension e) {
            throw new RuntimeException(e.getMessage());
        }
    }

    public double logPdf(double[] x) {
        return logPdf(x, mean, precision, logDet, 1.0);
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

        return dim * logNormalize + 0.5 * (logDet - dim*Math.log(scale) - SSE / scale);   // There was an error here.
        // Variance = (scale * Precision^{-1})

    }

    /* Equal precision, independent dimensions */
//    public static final double logPdf(double[] x, double[] mean, double precision, double scale) {
//
//        final int dim = x.length;
//
//        double SSE = 0;
//        for (int i = 0; i < dim; i++) {
//            double delta = x[i] - mean[i];
//            SSE += delta * delta;
//        }
//
//        return dim * logNormalize + 0.5 * (dim * Math.log(precision) - Math.log(scale) - SSE * precision / scale);
//        // TODO May not be correct
//    }

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

    public static void main(String[] args) {
        double[] start = new double[] {1, 2};
        double[] stop  = new double[] {0, 0};
        double[][] precision = new double[][] { {2,0.5},{0.5,1} };
        double scale = 0.2;
        System.err.println("logPDF = "+ logPdf(start, stop, precision, Math.log(calculatePrecisionMatrixDeterminate(precision)), scale));

    }

    public static final double logNormalize = -0.5 * Math.log(2.0 * Math.PI);

}
