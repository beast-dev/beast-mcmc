package test.dr.math;

import dr.math.EmpiricalBayesPoissonSmoother;
import dr.math.MathUtils;
import dr.stats.DiscreteStatistics;

/**
 * @author Marc A. Suchard
 */
public class EmpiricalBayesPoissonSmootherTest extends MathTestCase {

    public void testSmootherMean() {

        MathUtils.setSeed(666);
        int length = 100;

        double[] in = new double[length];
        for (int i = 0; i < length; ++i) {
            in[i] = MathUtils.nextDouble() * 100.0;
        }

        doSmoothing(in);
    }

    public void testUnstableMean() {
        double[] in = new double[] {0, 1};
        doSmoothing(in);
    }

    private void doSmoothing(double[] in) {

        double meanX = DiscreteStatistics.mean(in);
        double varX = DiscreteStatistics.variance(in);
        System.err.println("Original mean: " + meanX);
        System.err.println("Original var : " + varX + "\n");

        double[] out = EmpiricalBayesPoissonSmoother.smooth(in);

        double meanY = DiscreteStatistics.mean(out);
        double varY = DiscreteStatistics.variance(out);
        System.err.println("Smoothed mean: " + meanY);
        System.err.println("Smoothed var : " + varY);
        System.err.println("");

        assertEquals(meanX, meanY, tolerance);
        assertTrue(varY <= varX);
    }

    private static final double tolerance = 10E-6;


}
