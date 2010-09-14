package test.dr.distibutions;

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

            double x = Math.random() * 10;

            invGaussian.setMean(M);
            invGaussian.setShape(S);

            //double pdf = 1.0 / (x * S * Math.sqrt(2 * Math.PI)) * Math.exp(-Math.pow(Math.log(x) - M, 2) / (2 * S * S));
            double pdf = Math.sqrt(S/(2.0 * Math.PI * x * x * x)) * Math.exp((-1.0 * S * Math.pow((x - M), 2))/(2 * M * M * x));
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

            assertEquals(M, invGaussian.mean(), 1e-10);
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

            assertEquals(S, invGaussian.getShape(), 1e-10);
        }
    }

    public void testCDFAndQuantile() {		
        invGaussian.setMean(1.0);
        invGaussian.setShape(351.7561121947152);
        double q = invGaussian.quantile(0.20811009197062338);
        assertEquals(0.20811009197062338, invGaussian.cdf(q), 3.0e-3);

        for (int i = 0; i < 10000; i++) {
            double M = 1.0;
            double S = Math.random() * 1000.0 + 0.01;
            invGaussian.setMean(M);
            invGaussian.setShape(S);

            double p = Math.random()*0.98 + 0.01;
            double quantile = invGaussian.quantile(p);
            //System.out.println(quantile + "\t" + p + "\t" + M + "\t" + S);

            double cdf = invGaussian.cdf(quantile);
            if(((int)S)==351) {
                assertEquals(p, cdf, 1.0e-2);
            }
            else {
                assertEquals(p, cdf, 1.0e-3);
            }
        }

        /* Test with examples using R */
        invGaussian.setMean(5);
        invGaussian.setShape(0.5);
        assertEquals(0.75, invGaussian.cdf(3.022232), 1e-5);

        invGaussian.setMean(1.0);
        invGaussian.setShape(17.418709855826197);
        double q2 =invGaussian.quantile(0.27959422055126726);
        double p_hat = invGaussian.cdf(q2);
        assertEquals(0.27959422055126726, p_hat, 1.0e-3);

        invGaussian.setMean(1.0);
        invGaussian.setShape(0.4078303443934461);
        assertEquals(0.05514379243099207, invGaussian.cdf(invGaussian.quantile(0.05514379243099207)), 1.0e-3);
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
            zeroFinder.setMaximumIterations(100);
            zeroFinder.evaluate();
            assertEquals(f.quantile(i), zeroFinder.getResult(), 1e-3);
        }
    }

    public void testCDFAndQuantile3() {
        double[] shapes = {0.010051836, 0.011108997, 0.01227734, 0.013568559, 0.014995577,
                0.016572675, 0.018315639, 0.020241911, 0.022370772, 0.024723526, 0.027323722,
                0.030197383, 0.03337327, 0.036883167, 0.040762204, 0.045049202, 0.049787068,
                0.05502322, 0.060810063, 0.067205513, 0.074273578, 0.082084999, 0.090717953,
                0.100258844, 0.110803158, 0.122456428, 0.135335283, 0.149568619, 0.165298888,
                0.182683524, 0.201896518, 0.22313016, 0.246596964, 0.272531793, 0.301194212,
                0.332871084, 0.367879441, 0.40656966, 0.449328964, 0.496585304, 0.548811636,
                0.60653066, 0.670320046, 0.740818221, 0.818730753, 0.904837418, 1, 1.105170918,
                1.221402758, 1.349858808, 1.491824698, 1.648721271, 1.8221188, 2.013752707,
                2.225540928, 2.459603111, 2.718281828, 3.004166024, 3.320116923, 3.669296668,
                4.055199967, 4.48168907, 4.953032424, 5.473947392, 6.049647464, 6.685894442,
                7.389056099, 8.166169913, 9.025013499, 9.974182455, 11.02317638, 12.18249396,
                13.46373804, 14.87973172, 16.44464677, 18.17414537, 20.08553692, 22.19795128,
                24.5325302, 27.11263892, 29.96410005, 33.11545196, 36.59823444, 40.44730436,
                44.70118449, 49.40244911, 54.59815003, 60.3402876, 66.68633104, 73.6997937,
                81.45086866, 90.0171313, 99.48431564, 109.9471725, 121.5104175, 134.2897797,
                148.4131591, 164.0219073, 181.2722419, 200.33681, 221.4064162, 244.6919323,
                270.4264074, 298.867401, 330.2995599, 365.0374679, 403.4287935, 445.8577701,
                492.7490411, 544.5719101, 601.8450379, 665.141633, 735.0951892, 812.4058252,
                897.8472917, 992.2747156};
        for (double shape : shapes) {
            invGaussian.setShape(shape);
            for (double p = 0.01; p < 0.99; p += 0.01) {
                double q = invGaussian.quantile(p);
                double p_hat = invGaussian.cdf(q);
                assertEquals(p, p_hat, 1.0e-3);
            }
        }
    }
}
