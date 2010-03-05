package dr.math.distributions;

import dr.math.GammaFunction;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;

/**
 * @author Marc Suchard
 */
public class WishartDistribution implements MultivariateDistribution {

    public static final String TYPE = "Wishart";

    private int df;
    private int dim;
    private double[][] scaleMatrix;
    private Matrix S;
    private Matrix Sinv;
    private double logNormalizationConstant;

    /**
     * A Wishart distribution class for \nu degrees of freedom and scale matrix S
     * Expectation = \nu * S
     *
     * @param df
     * @param scaleMatrix
     */

    public WishartDistribution(int df, double[][] scaleMatrix) {
        this.df = df;
        this.scaleMatrix = scaleMatrix;
        this.dim = scaleMatrix.length;

        S = new Matrix(scaleMatrix);
        Sinv = S.inverse();

        computeNormalizationConstant();

    }

    public WishartDistribution(int dim) { // returns a non-informative (unormalizable) density
        this.df = 0;
        this.scaleMatrix = null;
        this.dim = dim;
        logNormalizationConstant = 1.0;
    }


    private void computeNormalizationConstant() {
        logNormalizationConstant = computeNormalizationConstant(new Matrix(scaleMatrix), df, dim);
    }

    public static double computeNormalizationConstant(Matrix Sinv, int df, int dim) {

        double logNormalizationConstant = 0;
        try {
            logNormalizationConstant = -df / 2.0 * Math.log(Sinv.determinant());
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }
        logNormalizationConstant -= df * dim / 2.0 * Math.log(2);
        logNormalizationConstant -= dim * (dim - 1) / 4.0 * Math.log(Math.PI);
        for (int i = 1; i <= dim; i++) {
            logNormalizationConstant -= GammaFunction.lnGamma((df + 1 - i) / 2.0);
        }
        return logNormalizationConstant;
    }


    public String getType() {
        return TYPE;
    }

    public double[][] getScaleMatrix() {
        return scaleMatrix;
    }

    public double[] getMean() {
        return null;
    }

    private void testMe() {

        int length = 100000;

        double save1 = 0;
        double save2 = 0;
        double save3 = 0;
        double save4 = 0;

        for (int i = 0; i < length; i++) {

            double[][] draw = nextWishart();
            save1 += draw[0][0];
            save2 += draw[0][1];
            save3 += draw[1][0];
            save4 += draw[1][1];

        }

        save1 /= length;
        save2 /= length;
        save3 /= length;
        save4 /= length;

        System.err.println("S1: " + save1);
        System.err.println("S2: " + save2);
        System.err.println("S3: " + save3);
        System.err.println("S4: " + save4);


    }

    public int df() {
        return df;
    }

    public double[][] scaleMatrix() {
        return scaleMatrix;
    }

    public double[][] nextWishart() {
        return nextWishart(df, scaleMatrix);
    }

    /**
     * Generate a random draw from a Wishart distribution
     * Follows Odell and Feiveson (1996) JASA 61, 199-203
     *
     * Returns a random variable with expectation = df * scaleMatrix
     *
     * @param df
     * @param scaleMatrix
     * @return
     */
    public static double[][] nextWishart(double df, double[][] scaleMatrix) {

        int dim = scaleMatrix.length;
        double[][] draw = new double[dim][dim];

        double[][] z = new double[dim][dim];

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < i; j++) {
                z[i][j] = MathUtils.nextGaussian();
            }
        }

        for (int i = 0; i < dim; i++)
            z[i][i] = Math.sqrt(MathUtils.nextGamma((df - i) * 0.5, 0.5));   // sqrt of chisq with df-i dfs

        double[][] cholesky = new double[dim][dim];
        for (int i = 0; i < dim; i++) {
            for (int j = i; j < dim; j++)
                cholesky[i][j] = cholesky[j][i] = scaleMatrix[i][j];
        }

        try {
            cholesky = (new CholeskyDecomposition(cholesky)).getL();
            // caution: this returns the lower triangular form
        } catch (IllegalDimension illegalDimension) {
        }

        double[][] result = new double[dim][dim];

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {     // lower triangular
                for (int k = 0; k < dim; k++)     // can also be shortened
                    result[i][j] += cholesky[i][k] * z[k][j];
            }
        }

        for (int i = 0; i < dim; i++) {           // lower triangular, so more efficiency is possible
            for (int j = 0; j < dim; j++) {
                for (int k = 0; k < dim; k++)
                    draw[i][j] += result[i][k] * result[j][k];   // transpose of 2nd element
            }
        }

        return draw;
    }

    public double logPdf(double[] x) {
        Matrix W = new Matrix(x, dim, dim);
        return logPdf(W, Sinv, df, dim, logNormalizationConstant);
    }

    public static double logPdf(Matrix W, Matrix Sinv, int df, int dim, double logNormalizationConstant) {

        double logDensity = 0;

//		System.err.println("yoyo "+df+" "+dim);

//        double det = 0;
        try {
            if (!W.isPD())
                return Double.NEGATIVE_INFINITY;

            final double det = W.determinant();
//		} catch (IllegalDimension illegalDimension) {
//			illegalDimension.printStackTrace();
//		}

            if (det <= 0) {
//			System.err.println("Not PD:\n"+W);
//			System.err.println("det = "+det);
//			System.exit(-1);
                return Double.NEGATIVE_INFINITY;
            }

//		try {
//			logDensity = Math.log(W.determinant());
            logDensity = Math.log(det);

            logDensity *= 0.5;
            logDensity *= df - dim - 1;

            // need only diagonal, no? seems a waste to compute
            // the whole matrix
            if (Sinv != null) {
                Matrix product = Sinv.product(W);

                for (int i = 0; i < dim; i++)
                    logDensity -= 0.5 * product.component(i, i);
            }

            // System.err.println(Sinv);

        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }

        logDensity += logNormalizationConstant;
        return logDensity;
    }

    public static void main(String[] argv) {
        WishartDistribution wd = new WishartDistribution(2, new double[][]{{500.0}});
        // The above is just an approximation
        GammaDistribution gd = new GammaDistribution(1.0 / 1000.0, 1000.0);
        double[] x = new double[]{1.0};
        System.out.println("Wishart, df=2, scale = 500, PDF(1.0): " + wd.logPdf(x));
        System.out.println("Gamma, shape = 1/1000, scale = 1000, PDF(1.0): " + gd.logPdf(x[0]));

        wd = new WishartDistribution(4, new double[][]{{5.0}});
        gd = new GammaDistribution(2.0, 10.0);
        x = new double[]{1.0};
        System.out.println("Wishart, df=4, scale = 5, PDF(1.0): " + wd.logPdf(x));
        System.out.println("Gamma, shape = 1/1000, scale = 10, PDF(1.0): " + gd.logPdf(x[0]));
        // These tests show the correspondence between a 1D Wishart and a Gamma

        wd = new WishartDistribution(1);
        x = new double[]{0.1};
        System.out.println("Wishart, uninformative, PDF(0.1): " + wd.logPdf(x));
        x = new double[]{1.0};
        System.out.println("Wishart, uninformative, PDF(1.0): " + wd.logPdf(x));
        x = new double[]{10.0};
        System.out.println("Wishart, uninformative, PDF(10.0): " + wd.logPdf(x));
        // These tests show the correspondence between a 1D Wishart and a Gamma
    }
}
