package test.dr.inference.timeseries;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Aggregates all time-series unit tests into a single runnable suite.
 * <p>
 * Run from IntelliJ or via:
 * <pre>
 *   java -cp build/beast.jar junit.textui.TestRunner test.dr.inference.timeseries.AllTimeSeriesTests
 * </pre>
 */
public class AllTimeSeriesTests {

    public static Test suite() {
        TestSuite suite = new TestSuite("Time-series module tests");
        suite.addTest(UniformTimeGridTest.suite());
        suite.addTest(OUProcessModelTest.suite());
        suite.addTest(GaussianObservationModelTest.suite());
        suite.addTest(KalmanLikelihoodEngineTest.suite());
        suite.addTest(KalmanGradientEngineTest.suite());
        suite.addTest(AnalyticalKalmanGradientEngineTest.suite());
        return suite;
    }

    public static void main(String[] args) {
        junit.textui.TestRunner.run(suite());
    }
}