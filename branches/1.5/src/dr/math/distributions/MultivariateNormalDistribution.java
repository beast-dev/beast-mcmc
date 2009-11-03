package dr.math.distributions;

import dr.math.MathUtils;
import dr.math.matrixAlgebra.*;

/**
 * @author Marc Suchard
 */
public class MultivariateNormalDistribution implements MultivariateDistribution {

    public static final String TYPE = "MultivariateNormal";

    private final double[] mean;
    private final double[][] precision;
    //	private int dim;
    private final double logDet;

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

    public static double calculatePrecisionMatrixDeterminate(double[][] precision) {
        try {
            return new Matrix(precision).determinant();
        } catch (IllegalDimension e) {
            throw new RuntimeException(e.getMessage());
        }
    }

    public double logPdf(double[] x) {
        return logPdf(x, mean, precision, logDet, 1.0);
    }

    public static double logPdf(double[] x, double[] mean, double[][] precision,
                                      double logDet, double scale) {

        if (logDet == Double.NEGATIVE_INFINITY)
            return logDet;

        final int dim = x.length;
        final double[] delta = new double[dim];
        final double[] tmp = new double[dim];

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
    public static double logPdf(double[] x, double[] mean, double precision, double scale) {

        final int dim = x.length;

        double SSE = 0;
        for (int i = 0; i < dim; i++) {
            double delta = x[i] - mean[i];
            SSE += delta * delta;
        }

        return dim * logNormalize + 0.5 * (dim * (Math.log(precision) - Math.log(scale)) - SSE * precision / scale);
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

        double[][] cholesky;

        try {
            cholesky = (new CholeskyDecomposition(variance)).getL();
        } catch (IllegalDimension illegalDimension) {
            throw new RuntimeException("Attempted Cholesky decomposition on non-square matrix");
        }

        double[] x = new double[dim];
        System.arraycopy(mean,0,x,0,dim);

        double[] epsilon = new double[dim];
        for (int i = 0; i < dim; i++)
            epsilon[i] = MathUtils.nextGaussian();

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j <= i; j++) {
                x[i] += cholesky[i][j] * epsilon[j];
                // caution: decomposition returns lower triangular
            }
        }
        return x;
    }

    // todo should be a test, no?
    public static void main(String[] args) {
        double[] start = {1, 2};
        double[] stop  = {0, 0};
        double[][] precision = { {2,0.5},{0.5,1} };
        double scale = 0.2;
        System.err.println("logPDF = "+ logPdf(start, stop, precision, Math.log(calculatePrecisionMatrixDeterminate(precision)), scale));
        System.err.println("Should = -19.94863\n");

        System.err.println("logPDF = "+logPdf(start,stop,2,0.2));
        System.err.println("Should = -24.53529\n");

        System.err.println("Random draws: ");
        int length = 10000;
        double[] mean = new double[2];
        double[] SS   = new double[2];
        double[] var  = new double[2];
        double ZZ = 0;
        for(int i=0; i<length; i++) {
            double[] draw = nextMultivariateNormalPrecision(start,precision);
            for(int j=0; j<2; j++) {
                mean[j] += draw[j];
                SS[j] += draw[j]*draw[j];
            }
            ZZ += draw[0]*draw[1];
        }

        for(int j=0; j<2; j++) {
            mean[j] /= length;
            SS[j] /= length;
            var[j] = SS[j] - mean[j]*mean[j];
        }
        ZZ /= length;
        ZZ -= mean[0]*mean[1];

        System.err.println("Mean: "+new Vector(mean));
        System.err.println("TRUE: [ 1 2 ]\n");
        System.err.println("MVar: "+new Vector(var));
        System.err.println("TRUE: [ 0.571 1.14 ]\n");
        System.err.println("Covv: "+ZZ);
        System.err.println("TRUE: -0.286");

    }

    public static final double logNormalize = -0.5 * Math.log(2.0 * Math.PI);

}
