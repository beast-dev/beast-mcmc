package dr.math.distributions;

import dr.math.GammaFunction;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;

/**
 * @author Marc Suchard
 */
public class InverseWishartDistribution implements MultivariateDistribution {

    public static final String TYPE = "InverseWishart";

    private int df;
    private int dim;
    private double[][] scaleMatrix;
    private Matrix S;
    private double logNormalizationConstant;

    /**
     * An Inverser Wishart distribution class for \nu degrees of freedom and scale matrix S with dim k
     * Expectation = (\nu - k - 1)^{-1} * S
     *
     * @param df
     * @param scaleMatrix
     */

    public InverseWishartDistribution(int df, double[][] scaleMatrix) {
        this.df = df;
        this.scaleMatrix = scaleMatrix;
        this.dim = scaleMatrix.length;

        S = new Matrix(scaleMatrix);
        computeNormalizationConstant();
    }

    private void computeNormalizationConstant() {
        logNormalizationConstant = 0;
        try {
            logNormalizationConstant = df / 2.0 * Math.log(new Matrix(scaleMatrix).determinant());
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }
        logNormalizationConstant -= df * dim / 2.0 * Math.log(2);
        logNormalizationConstant -= dim * (dim - 1) / 4.0 * Math.log(Math.PI);
        for (int i = 1; i <= dim; i++) {
            logNormalizationConstant -= GammaFunction.lnGamma((df + 1 - i) / 2.0);
        }

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

    public int df() {
        return df;
    }

    public double[][] scaleMatrix() {
        return scaleMatrix;
    }


    public double logPdf(double[] x) {
        Matrix W = new Matrix(x, dim, dim);
        double logDensity = 0;

//	    System.err.println("here");
//	    double det = 0;
//	    try {
//	        det = W.determinant();
//	    }   catch (IllegalDimension illegalDimension) {
//            illegalDimension.printStackTrace();
//        }
//	    if( det < 0 ) {
//		    System.err.println("not positive definite");
//		    return Double.NEGATIVE_INFINITY;
//	    }


        try {
            logDensity = Math.log(W.determinant());

            logDensity *= -0.5;
            logDensity *= df + dim + 1;

            Matrix product = S.product(W.inverse());

            for (int i = 0; i < dim; i++)
                logDensity -= 0.5 * product.component(i, i);

        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }

        logDensity += logNormalizationConstant;
        return logDensity;
    }

}
