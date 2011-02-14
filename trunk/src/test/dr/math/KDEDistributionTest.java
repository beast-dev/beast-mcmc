package test.dr.math;

import dr.math.MathUtils;
import dr.math.distributions.NormalKDEDistribution;
import dr.math.matrixAlgebra.Vector;

/**
 * @author Marc A. Suchard
 */
public class KDEDistributionTest extends MathTestCase {

    public void testNormalKDE() {

        int length = 100;
        Double[] sample = new Double[length];
        for (int i = 0; i < length; i++) {
            sample[i] = (double) (i + 1);
        }

//        final int gridSize = 256;

        NormalKDEDistribution kde = new NormalKDEDistribution(sample, null, null, null);
        System.out.println("bw = " + kde.getBandWidth());
        assertEquals(rBandWidth[0], kde.getBandWidth(), tolerance);


        final int gridSize = NormalKDEDistribution.MINIMUM_GRID_SIZE;
        double[] testPoints = new double[gridSize];
        double x = kde.getFromPoint();
        double delta = (kde.getToPoint() - kde.getFromPoint()) / (gridSize - 1);
        for (int i = 0; i < gridSize; ++i) {
            testPoints[i] = x;
            x += delta;
        }
        System.err.println("Eval @ " + new Vector(testPoints));

        double[] testDensity = new double[gridSize];
        for (int i = 0; i < gridSize; ++i) {
            testDensity[i] = kde.pdf(testPoints[i]);
        }
        System.err.println("Den    " + new Vector(testDensity));

        System.err.println("den[0] = " + testDensity[0]);
        System.err.println("den[N] = " + testDensity[NormalKDEDistribution.MINIMUM_GRID_SIZE - 1]);

      //  System.exit(-1);
    }

    public void testGammaKDE() {

        int length = 10000;
        double shape = 4;
        double scale = 5;

        double[] values = new double[length];
        for (int i = 0; i < length; i++) {
            values[i] = MathUtils.nextGamma(shape, scale);
        }


//        GammaKDEDistribution kde = new GammaKDEDistribution(values);
//             System.err.println("prediction: at 2.02: "+kde.pdf(2.02177)+" at 0.405: "+kde.pdf(0.4046729)+" at 0.15: "+kde.pdf(0.1502078));
//     System.err.println("sm: "+kde.sampleMean());
        // TODO Need test values

    }

    private static double[] rBandWidth = { 12.24266 };
    private static double tolerance = 1E-4;
}

/*

# R test

sample = seq(1:100)
d = density(sample, bw="nrd")


*/