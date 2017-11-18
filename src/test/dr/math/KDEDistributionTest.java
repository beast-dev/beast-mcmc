package test.dr.math;

import dr.math.MathUtils;
import dr.math.distributions.NormalKDEDistribution;
import dr.math.matrixAlgebra.Vector;

/**
 * @author Marc A. Suchard
 */
public class KDEDistributionTest extends MathTestCase {
    
    public void testNormalCDF() {

        MathUtils.setSeed(666);
        int length = 100000;

        Double[] values = new Double[length];
        for (int i = 0; i < length; ++i) {
            values[i] = MathUtils.nextGaussian();
        }

        NormalKDEDistribution distribution = new NormalKDEDistribution(values);

        double cdf = distribution.cdf(0);
        System.out.println("cdf(0.00) = " + cdf);
        assertEquals(cdf, 0.5, cdfTolerance);

        cdf = distribution.cdf(1.96);
        System.out.println("cdf(1.96) = " + cdf);
        assertEquals(cdf, 0.975, cdfTolerance);

        double quantile = distribution.quantile(0.68);
        System.out.println("quantile(0.68) = " + quantile);
        assertEquals(quantile, 0.468, cdfTolerance);

        quantile = distribution.quantile(0.25);
        System.out.println("quantile(0.25) = " + quantile);
        assertEquals(quantile, -0.674, cdfTolerance);
    }

    private final static double cdfTolerance = 1E-2;

    public void testNormalKDE() {

        int length = 100;
        Double[] sample = new Double[length];
        for (int i = 0; i < length; i++) {
            sample[i] = (double) (i + 1);
        }

//        final int gridSize = 256;

        NormalKDEDistribution kde = new NormalKDEDistribution(sample, null, null, null);
        System.out.println("bw = " + kde.getBandWidth());
        double tolerance = 1E-4;
        assertEquals(rBandWidth[0], kde.getBandWidth(), tolerance);


        final int gridSize = NormalKDEDistribution.MINIMUM_GRID_SIZE;
        double[] testPoints = new double[gridSize];
        double x = kde.getFromPoint();
        double delta = (kde.getToPoint() - kde.getFromPoint()) / (gridSize - 1);
        for (int i = 0; i < gridSize; ++i) {
            testPoints[i] = x;
            x += delta;
        }
        System.out.println("Eval @ " + new Vector(testPoints));

        double[] testDensity = new double[gridSize];
        for (int i = 0; i < gridSize; ++i) {
            testDensity[i] = kde.pdf(testPoints[i]);
        }
        System.out.println("Den    " + new Vector(testDensity));

        System.out.println("den[0] = " + testDensity[0]);
        System.out.println("den[N] = " + testDensity[NormalKDEDistribution.MINIMUM_GRID_SIZE - 1]);

      //  System.exit(-1);
    }

    public void testGammaKDE() {

//        int length = 10000;
//        double shape = 4;
//        double scale = 5;

//        double[] values = new double[length];
//        for (int i = 0; i < length; i++) {
//            values[i] = MathUtils.nextGamma(shape, scale);
//        }


//        GammaKDEDistribution kde = new GammaKDEDistribution(values);
//             System.out.println("prediction: at 2.02: "+kde.pdf(2.02177)+" at 0.405: "+kde.pdf(0.4046729)+" at 0.15: "+kde.pdf(0.1502078));
//     System.out.println("sm: "+kde.sampleMean());
        // TODO Need test values

    }

    private static double[] rBandWidth = { 12.24266 };
}

/*

# R test

sample = seq(1:100)
d = density(sample, bw="nrd")


*/