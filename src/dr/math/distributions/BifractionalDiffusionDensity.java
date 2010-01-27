package dr.math.distributions;

import dr.math.UnivariateFunction;
import cern.jet.stat.Gamma;


/**
 * @author Marc Suchard
 *         <p/>
 *         **
 *         Most of the following functions come from:
 *         <p/>
 *         R Gorenflo, A Iskenderov and Y Luchko (2000) Mapping between solutions of fractional diffusion-wave equations.
 *         Fractional Calculus and Applied Analysis, 3, 75 - 86
 *         <p/>
 *         Also see:
 *         F. Mainardi and G. Pagnini (2003) The Wright functions as solutions of the time-fractional diffusion equation.
 *         Applied Mathematics and Computation, 141, pages?
 */

public class BifractionalDiffusionDensity implements Distribution {

    public BifractionalDiffusionDensity(double alpha, double beta) {
        this.alpha = alpha;
        this.beta = beta;
    }

    /**
     * probability density function of the distribution
     *
     * @param x argument
     * @return pdf value
     */
    public double pdf(double x) {
        return Math.exp(logPdf(x));
    }

    /**
     * the natural log of the probability density function of the distribution
     *
     * @param x argument
     * @return log pdf value
     */
    public double logPdf(double x) {
        return logPdf(x, 1.0, alpha, beta);
    }

    /**
     * cumulative density function of the distribution
     *
     * @param x argument
     * @return cdf value
     */
    public double cdf(double x) {
        throw new RuntimeException("Not yet implemented");
    }

    /**
     * quantile (inverse cumulative density function) of the distribution
     *
     * @param y argument
     * @return icdf value
     */
    public double quantile(double y) {
        throw new RuntimeException("Not yet implemented");
    }

    /**
     * mean of the distribution
     *
     * @return mean
     */
    public double mean() {
        throw new RuntimeException("Not yet implemented");
    }

    /**
     * variance of the distribution
     *
     * @return variance
     */
    public double variance() {
        throw new RuntimeException("Not yet implemented");
    }

    /**
     * @return a probability density function representing this distribution
     */
    public UnivariateFunction getProbabilityDensityFunction() {
        throw new RuntimeException("Not yet implemented");
    }

   /**
     * Taken from:  Saichev AI and Zaslavsky GM (1997) Fractional kinetic equations: solutions and applications.
     *              Chaos, 7, 753-764
     * @param x evaluation point
     * @param t evaluation time
     * @param alpha coefficient
     * @param beta coefficient
     * @return probability density
     */

    public static double logPdf(double x, double t, double alpha, double beta) {
        final double mu = beta / alpha;
        final double absX = Math.abs(x);
        final double absY = absX / Math.pow(t,mu);
        double density = 0;
        double incr = Double.MAX_VALUE;
        int m = 1; // 0th term = 0 \propto cos(pi/2)
        int sign = -1;

        while (Math.abs(incr) > eps && m < maxK) {
            incr =  sign / Math.pow(absY, m * alpha)
                         * Gamma.gamma(m * alpha + 1)
                         / Gamma.gamma(m * beta + 1)
                         * Math.cos(halfPI * (m * alpha + 1));
            density += incr;
            sign *= -1;
            m++;
        }

        return Math.log(density / (Math.PI * absX));
    }


    static class SignedDouble {
        double x;
        boolean positive;

        SignedDouble(double x, boolean positive) {
            this.x = x;
            this.positive = positive;
        }
    }

    private static SignedDouble logGamma(double z) {

        // To extend the gamma function to negative (non-integer) numbers, apply the relationship
        // \Gamma(z) = \frac{ \Gamma(z+n) }{ z(z+1)\cdots(z+n-1),
        // by choosing n such that z+n is positive
        if (z > 0) {
            return new SignedDouble(Gamma.logGamma(z), true);
        }
        int n = ((int) -z) + 1;
        boolean positive = (n % 2 == 0);
        return new SignedDouble(Gamma.logGamma(z + n) - Gamma.logGamma(-z + 1) + Gamma.logGamma(-z - n + 1), positive);
    }


    private static double gamma(double z) {

        // To extend the gamma function to negative (non-integer) numbers, apply the relationship
        // \Gamma(z) = \frac{ \Gamma(z+n) }{ z(z+1)\cdots(z+n-1),
        // by choosing n such that z+n is positive
        if (z > 0.0) {
            return Gamma.gamma(z);
        }
        int n = ((int) -z);
        if (z+n == 0.0) {
            return Double.NaN;
        }
        n++;
        boolean positive = (n % 2 == 0);
        double result = Gamma.gamma(z + n) / Gamma.gamma(-z + 1) * Gamma.gamma(-z - n + 1);
        if (!positive) result *= -1;
        return result;
    }

    public static final int maxK = 50;

    public static double infiniteSumBetaGreaterThanAlpha1(double z, double alpha, double beta) {
        double sum = 0;
        int k = 0;
        boolean isPositive = true;

        double incr = Double.MAX_VALUE;
        while (//(Math.abs(incr) > 1E-20) &&
                (k < maxK)) {

            double x1, x2, x3;
            x1 = gamma(0.5 - 0.5 * alpha - 0.5 * alpha * k);
            if (!Double.isNaN(x1)) {
                incr = x1;
            } else {
                System.err.println("Big problem!");
                System.exit(-1);
            }
            x2 = gamma(1.0 - beta - beta * k);
            if (!Double.isNaN(x2)) {
                incr /= x2;
            } else {
                 incr = 0;
            }
            x3 = gamma(0.5 * alpha + 0.5 * alpha * k);
            if (!Double.isNaN(x3)) {
                incr /= x3;
            } else {
                incr = 0;
            }
            incr *= Math.pow(z, k);
            if (isPositive) {
                sum += incr;
            } else {
                sum -= incr;
            }
            isPositive = !isPositive;
            k++;
        }
//        System.err.println("sum 1 Steps = "+k);
//        System.err.println("Last incr = "+incr);
        return sum;
    }

    public static double infiniteSumBetaGreaterThanAlpha2(double z, double alpha, double beta) {
        double sum = 0;
        int m = 0;
        boolean isPositive = true;

        double incr = Double.MAX_VALUE;
        while (// (Math.abs(incr) > 1E-20) &&
                (m < maxK)) {

            double x1, x2, x3, x4;

            x1 = gamma(1.0 / alpha + 2.0 / alpha * m);
            if (!Double.isNaN(x1)) {
                incr = x1;
            } else {
                System.err.println("Big problem!");
                System.exit(-1);
            }
            x2 = gamma(1.0 - 1.0 / alpha - 2.0 / alpha * m);
            if (!Double.isNaN(x2)) {
                incr *= x2;
            } else {
                System.err.println("Big problem!");
                System.exit(-1);
            }
            x3 = gamma(0.5 + m);
            if (!Double.isNaN(x3)) {
                incr /= x3;
            } else {
                incr = 0;
            }
            x4 = gamma(1.0 - beta / alpha - 2 * beta / alpha * m);
            if (!Double.isNaN(x4)) {
                incr /= x4;
            } else {
                incr = 0;
            }
            incr /= gamma(m + 1);
            incr *= Math.pow(z, m);
            if (isPositive) {
                sum += incr;
            } else {
                sum -= incr;
            }
            isPositive = !isPositive;
            m++;

        }
//        System.err.println("sum2 Steps = " + m);
//        System.err.println("Last incr = " + incr);
        return sum;
    }

    public static double evaluateGreensFunctionAlphaGreaterThanBeta(double x, double t, double alpha, double beta) {
        double z = Math.pow(2.0, alpha) * Math.pow(t, beta) / Math.pow(Math.abs(x), alpha);

        double[][] aAp = new double[][] { {0.5, alpha / 2.0}, {1.0, 1.0} };
        double[][] bBp = new double[][] { {1.0, beta}, {0.0, -alpha / 2.0} };

        return 1.0 / (Math.sqrt(Math.PI) * Math.abs(x)) * generalizedWrightFunction(-z, aAp, bBp);
    }


    public static double evaluateGreensFunctionBetaGreaterThanAlpha(double x, double t, double alpha, double beta) {

        double z1 = Math.pow(Math.abs(x),alpha) / (Math.pow(2,alpha) * Math.pow(t, beta));
        double z2 = Math.pow(x, 2.0) / (4.0 * Math.pow(t, 2 * beta / alpha));

        // Using specialized functions
        double green1 = 1.0 / Math.sqrt(Math.PI) * Math.pow(Math.abs(x), alpha - 1.0) /
                (Math.pow(2,alpha) * Math.pow(t, beta)) *
                infiniteSumBetaGreaterThanAlpha1(z1, alpha, beta);
        double green2 = 1.0 / Math.sqrt(Math.PI) /  (alpha * Math.pow(t, beta/alpha)) *
                infiniteSumBetaGreaterThanAlpha2(z2, alpha, beta);

        // Using a general function
        double[][] aAp1 = new double[][] { {0.5 - alpha / 2.0, -alpha / 2.0}, {1.0, 1.0} };
        double[][] bBp1 = new double[][] { {1.0 - beta, -beta}, {alpha / 2.0, alpha / 2.0} };

        double[][] aAp2 = new double[][] { {1.0 / alpha, 2.0 / alpha}, {1.0 - 1.0 / alpha, -2.0 / alpha} };
        double[][] bBp2 = new double[][] { {0.5, 1.0 }, {1.0 - beta / alpha, -2.0 * beta / alpha} };

        double green3 = 1.0 / Math.sqrt(Math.PI) * Math.pow(Math.abs(x), alpha - 1.0) /
                (Math.pow(2,alpha) * Math.pow(t, beta)) *
                generalizedWrightFunction(-z1, aAp1, bBp1);
        double green4 = 1.0 / Math.sqrt(Math.PI) /  (alpha * Math.pow(t, beta/alpha)) *
                generalizedWrightFunction(-z2, aAp2, bBp2);

        // Compare two methods
        if (Math.abs(green1 + green2 - green3 - green4) > 1E-10) {
            System.err.println("Computation error");
            System.exit(-1);
        }

        return green1 + green2;
    }

    public static double generalizedWrightFunction(double z, double[][] aAp, double[][] bBq) {
        final int p = aAp.length;
        final int q = bBq.length;
        double sum = 0;
        double incr;
        int k = 0;

        while (// incr > eps &&
                k < maxK) {
            incr = 1;
            for (int i = 0; i < p; i++) {
                final double[] aAi = aAp[i];

                double x = gamma(aAi[0] + aAi[1] * k); // TODO Precompute these factors
                if (!Double.isNaN(x)) {
                    incr *= x;
                } else {
                    incr = Double.POSITIVE_INFINITY;
                }
            }
            for (int i = 0; i < q; i++) {
                final double[] bBi = bBq[i];
                double x = gamma(bBi[0] + bBi[1] * k); // TODO Precompute these factors
                if (!Double.isNaN(x)) {
                    incr /= x;
                } else {
                    incr = 0.0;
                }
            }
            incr /= gamma(k+1); // k! TODO Precompute these factors
            incr *= Math.pow(z,k); // TODO Just multiply out

            sum += incr;
            k++;
        }
        return sum;
    }

    public static void main(String[] arg) {

        double alpha = 2.0;
        double beta = 0.8;
        double z = 2.0;
        double z1 = -2.34;
        SignedDouble result = logGamma(z1);
        System.err.println("logGamma("+z1+") = "+result.x+" "+(result.positive ? "(+)" : "(-)"));
        System.err.println("gamma("+z1+") = "+ gamma(z1));
        System.err.println("gamma(-2.0) = "+gamma(-2.0));
        System.err.println("Result1 = "+infiniteSumBetaGreaterThanAlpha1(z, alpha, beta));
        System.err.println("Result2 = "+infiniteSumBetaGreaterThanAlpha2(z, alpha, beta));
        System.err.println("");

        double var = 4.0;
        double t = 0.5 * var;
        double x = 1.0;
        System.err.println("p(x = "+x+", v = "+var+") = "+ evaluateGreensFunctionBetaGreaterThanAlpha(x,t,alpha, beta));


        System.err.println("");
    }
    

    private static double eps = 1E-20;
    private static double halfPI = Math.PI / 2.0;

    private double alpha;
    private double beta;
}

/*

# R script

G = function(x = 0, t, alpha, beta) {
  density = NaN
  if (x == 0) {
  	if (beta == 1) {
  		density = gamma(1/alpha) / (pi * alpha * (t^(1/alpha)))
  	} else {
  		density = 1 / (alpha * t^(beta/alpha)) / sin(pi / alpha) / gamma(1 - (beta/alpha))
  	}
  }
  density
}


*/