package dr.math.distributions;

import dr.math.GammaFunction;

/**
 * @author Marc A. Suchard
 */
public class DirichletDistribution implements MultivariateDistribution {

    public static final String TYPE = "dirichletDistribution";

    private double[] counts;
    private double countSum = 0.0;
    private int dim;

    private double logNormalizingConstant;

    public DirichletDistribution(double[] counts) {
        this.counts = counts;
        dim = counts.length;
        for (int i = 0; i < dim; i++)
            countSum += counts[i];

        computeNormalizingConstant();
    }

    private void computeNormalizingConstant() {
        logNormalizingConstant = GammaFunction.lnGamma(countSum);
        for (int i = 0; i < dim; i++)
            logNormalizingConstant -= GammaFunction.lnGamma(counts[i]);
    }


    public double logPdf(double[] x) {

        if (x.length != dim) {
            throw new IllegalArgumentException("data array is of the wrong dimension");
        }

        double logPDF = logNormalizingConstant;
        for (int i = 0; i < dim; i++) {
            logPDF += (counts[i] - 1) * Math.log(x[i]);
            if (x[i] <= 0.0 || x[i] >= 1.0) {
                logPDF = Double.NEGATIVE_INFINITY;
                break;
            }
        }       
        return logPDF;
    }

    public double[][] getScaleMatrix() {
        return null;
    }

    public double[] getMean() {
        double[] mean = new double[dim];
        for (int i = 0; i < dim; i++)
            mean[i] = counts[i] / countSum;
        return mean;
    }

    public String getType() {
        return TYPE;
    }
}
