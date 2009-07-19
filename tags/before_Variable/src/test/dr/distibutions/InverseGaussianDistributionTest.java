package test.dr.distibutions;

//import dr.math.distributions.LogNormalDistribution;
import dr.math.distributions.InverseGaussianDistribution;
import dr.math.interfaces.OneVariableFunction;
import dr.math.iterations.BisectionZeroFinder;
import junit.framework.TestCase;

/**
 * @author Wai Lok Sibon Li
 */
public class InverseGaussianDistributionTest extends TestCase {

    InverseGaussianDistribution invGaussian;

    public void setUp() {

        invGaussian = new InverseGaussianDistribution(1.0, 2.0);
    }

    public void testPdf() {

        System.out.println("Testing 10000 random pdf calls");

        for (int i = 0; i < 10000; i++) {
            double M = Math.random() * 10.0 + 0.1;
            double S = Math.random() * 5.0 + 0.01;

            //double x = Math.log(Math.random() * 10);
            double x = Math.random() * 10;

            invGaussian.setMean(M);
            invGaussian.setShape(S);

            //double pdf = 1.0 / (x * S * Math.sqrt(2 * Math.PI)) * Math.exp(-Math.pow(Math.log(x) - M, 2) / (2 * S * S));
            double pdf = Math.sqrt(S/(2.0 * Math.PI * x * x * x)) * Math.exp((-1.0 * S * Math.pow((x - M), 2))/(2 * M * M * x));

            //System.out.println("Testing invGaussian[M=" + M + " S=" + S + "].pdf(" + x + ")");

            assertEquals(pdf, invGaussian.pdf(x), 1e-10);
        }

        /* Test with an example using R */
        invGaussian.setMean(2.835202292812448);
        invGaussian.setShape(3.539139491639669);
        assertEquals(0.1839934, invGaussian.pdf(2.540111), 1e-6);
    }

    public void testMean() {

        for (int i = 0; i < 1000; i++) {
            double M = Math.random() * 10.0 + 0.1;
            double S = Math.random() * 5.0 + 0.01;

            invGaussian.setMean(M);
            invGaussian.setShape(S);

            double mean = M;

            //System.out.println("Testing invGaussian[M=" + M + " S=" + S + "].mean()");

            assertEquals(mean, invGaussian.mean(), 1e-10);
        }
    }

    public void testVariance() {

        for (int i = 0; i < 1000; i++) {
            double M = Math.random() * 10.0 + 0.1;
            double S = Math.random() * 5.0 + 0.01;

            invGaussian.setMean(M);
            invGaussian.setShape(S);
            double variance = (M * M * M) / S;
            assertEquals(variance, invGaussian.variance(), 1e-8);
        }
    }

    public void testShape() {

        System.out.println("Testing 10000 random quantile(0.5) calls");

        for (int i = 0; i < 10000; i++) {
            double M = Math.random() * 10.0 + 0.1;
            double S = Math.random() * 5.0 + 0.01;

            invGaussian.setMean(M);
            invGaussian.setShape(S);

            double shape = S;
            assertEquals(shape, invGaussian.getShape(), 1e-10);
        }
    }

    public void testCDFAndQuantile() {

        System.out.println("Testing 100000 random quantile/cdf pairs");

        for (int i = 0; i < 10000; i++) {

            double M = Math.random() * 10.0 + 0.5;
            double S = Math.random() * 5.0 + 0.01;

            invGaussian.setMean(M);
            invGaussian.setShape(S);

            double p = Math.random();
            double quantile = invGaussian.quantile(p);

            double cdf = invGaussian.cdf(quantile);
            //System.out.println("p (input) " + p+"\tquantile " + quantile + "\toutput " + cdf + "\t" + M + "\t" + S);
            assertEquals(p, cdf, 1e-5);
        }

        /* Test with an example using R */
        invGaussian.setMean(5);
        invGaussian.setShape(0.5);
        assertEquals(0.75, invGaussian.cdf(3.022232), 1e-5);
    }

    public void testCDFAndQuantile2() {

        final InverseGaussianDistribution f = new InverseGaussianDistribution(1, 1);
        for (double i = 0.01; i < 0.95; i += 0.01) {
            final double y = i;

            BisectionZeroFinder zeroFinder = new BisectionZeroFinder(new OneVariableFunction() {
                public double value(double x) {
                    return f.cdf(x) - y;
                }
            }, 0.01, 100);
            zeroFinder.evaluate();
            assertEquals(f.quantile(i), zeroFinder.getResult(), 1e-6);
        }
    }
}
