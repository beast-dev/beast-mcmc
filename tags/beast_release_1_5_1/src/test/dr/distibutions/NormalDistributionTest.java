package test.dr.distibutions;

import dr.math.distributions.NormalDistribution;
import dr.math.interfaces.OneVariableFunction;
import dr.math.iterations.BisectionZeroFinder;
import junit.framework.TestCase;

/**
 * @author Wai Lok Sibon Li
 * 
 */
public class NormalDistributionTest extends TestCase {
    NormalDistribution norm;

    public void setUp() {

        norm = new NormalDistribution(0.0, 1.0);
    }


    public void testPdf() {

        System.out.println("Testing 10000 random pdf calls");

        for (int i = 0; i < 10000; i++) {
            double M = Math.random() * 10.0 - 5.0;
            double S = Math.random() * 10;

            double x = Math.random() * 10;

            norm.setMean(M);
            norm.setSD(S);
            
            double a = 1.0 / (Math.sqrt(2.0 * Math.PI) * S);
            double b = -(x - M) * (x - M) / (2.0 * S * S);
            double pdf =  a * Math.exp(b);

            assertEquals(pdf, norm.pdf(x), 1e-10);
        }

        /* Test with an example using R */
        norm.setMean(2.835202292812448);
        norm.setSD(3.539139491639669);
        assertEquals(0.1123318, norm.pdf(2.540111), 1e-6);
    }

    public void testMean() {

        for (int i = 0; i < 1000; i++) {
            double M = Math.random() * 10.0 - 5.0;

            norm.setMean(M);

            assertEquals(M, norm.mean(), 1e-10);
        }
    }

    public void testVariance() {

        for (int i = 0; i < 1000; i++) {
            double S = Math.random() * 10;
            norm.setSD(S);

            double variance = S * S;

            assertEquals(variance, norm.variance(), 1e-10);
        }
    }


    public void testMedian() {

        System.out.println("Testing 10000 random quantile(0.5) calls");

        for (int i = 0; i < 10000; i++) {
            double M = Math.random() * 10.0 - 5.0;
            double S = Math.random() * 10;

            norm.setMean(M);
            norm.setSD(S);

            double median = M;

            assertEquals(median, norm.quantile(0.5), 1e6);
        }
    }

    public void testCDFAndQuantile() {

        System.out.println("Testing 10000 random quantile/cdf pairs");

        for (int i = 0; i < 10000; i++) {

            double M = Math.random() * 10.0 - 5.0;
            double S = Math.random() * 10;

            norm.setMean(M);
            norm.setSD(S);

            double p = Math.random();
            double quantile = norm.quantile(p);
            double cdf = norm.cdf(quantile);

            assertEquals(p, cdf, 1e-8);
        }

    }

    public void testCDFAndQuantile2() {
        for(int i=0; i<10000; i++) {
            double x =Math.random();
            double m = Math.random() * 10;
            double s = Math.random() * 10;
            
            double a = NormalDistribution.cdf(x, m, s, false);
            double b =NormalDistribution.cdf(x, m, s);
            
            assertEquals(a, b, 1.0e-8);
        }
    }


}
