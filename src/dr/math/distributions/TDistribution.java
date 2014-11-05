package dr.math.distributions;

import dr.math.GammaFunction;
import dr.math.UnivariateFunction;

/**
 * @author Marc A. Suchard
 */
public class TDistribution implements Distribution {

    public TDistribution(double center, double scale, double df) {
        this.center = center;
        this.scale = scale;
        this.df = df;
    }

    public double pdf(double x) {
        return Math.exp(logPDF(x, center, scale, df));
    }

    public double logPdf(double x) {
        return logPDF(x, center, scale, df);
    }

    public double cdf(double x) {
        throw new RuntimeException("Not yet implemented");
    }

    public double quantile(double y) {
        throw new RuntimeException("Not yet implemented");
    }

    public double mean() {
        if (df > 1)
            return center;
        return Double.NaN;
    }

    public double variance() {
        if (df > 2)
            return scale * df / (df - 2);
        return Double.NaN;
    }

    public UnivariateFunction getProbabilityDensityFunction() {
        throw new RuntimeException("Not yet implemented");
    }

    /*
          * t-distribution
          *
          *  \frac{ \Gamma(\frac{\nu+1}{2}) }{ \sqrt{\nu \pi} \Gamma(\frac{\nu}{2}) } *
          *      \left(
          *        1 + \frac{(x-x_0)^2}{\nu}
          *      \right)^{ -\frac{\nu+1}{2} }
          *
          * Variance = \frac{\nu}{\nu-2}
          *
          */

    /* Returns the log PDF of a t-distribution
      * @param df degrees of freedom (=1 cauchy, =\infty normal)
      * @param x0 mean for df > 1
      * @param scale variance = scale * df / (df - 2) for df > 2
      */

    public static double logPDF(double x, double x0, double scale, double df) {

        double loc = x - x0;

        double logPDF = GammaFunction.lnGamma((df + 1) / 2)
                - 0.5 * Math.log(df)
                - logSqrtPi
                - 0.25 * Math.log(scale)
                - GammaFunction.lnGamma(df / 2)
                - (df + 1) / 2 * Math.log(1 + loc * loc / df / scale);

        return logPDF;
    }

    private double df;
    private double scale;
    private double center;

    private static double logSqrtPi = Math.log(Math.sqrt(Math.PI));

}
