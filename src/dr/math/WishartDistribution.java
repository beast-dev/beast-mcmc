package dr.math;

import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.IllegalDimension;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Jun 12, 2007
 * Time: 9:01:16 PM
 * To change this template use File | Settings | File Templates.
 */
public class WishartDistribution implements MultivariateDistribution {

    public static final String TYPE = "Wishart";

    private int df;
    private double[][] inverseScaleMatrix;

    /**
     * A Wishart distribution class for \nu degrees of freedom and inverse scale matrix S
     * Expectation = \nu * S
     *
     * @param df
     * @param inverseScaleMatrix
     */

    public WishartDistribution(int df, double[][] inverseScaleMatrix) {
        this.df = df;
        this.inverseScaleMatrix = inverseScaleMatrix;
//		testMe();
    }

    public String getType() {
        return TYPE;
    }

    public double[][] getScaleMatrix() {
        return inverseScaleMatrix;
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

    public double[][] inverseScaleMatrix() {
        return inverseScaleMatrix;
    }

    public double[][] nextWishart() {
        return nextWishart(df, inverseScaleMatrix);
    }

    /**
     * Generate a random draw from a Wishart distribution
     * Follows Odell and Feiveson (1996) JASA 61, 199-203
     *
     * @param df
     * @param inverseScaleMatrix
     * @return
     */
    public static double[][] nextWishart(int df, double[][] inverseScaleMatrix) {

        int dim = inverseScaleMatrix.length;
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
                cholesky[i][j] = cholesky[j][i] = inverseScaleMatrix[i][j];
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
//		return result;
    }

    public double logPdf(Parameter x) {
//		MatrixParameter matrix = (MatrixParameter)x;
        // Currently Wishart distributions are only used in Gibbs steps, so calculating densities is not necessary
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
