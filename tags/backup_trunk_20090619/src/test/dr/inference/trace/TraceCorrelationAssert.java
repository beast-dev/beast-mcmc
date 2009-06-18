package test.dr.inference.trace;

import dr.inference.trace.TraceCorrelation;
import junit.framework.TestCase;

/**
 * @author Alexei Drummond
 */
public class TraceCorrelationAssert extends TestCase {

    public TraceCorrelationAssert(String name) {
        super(name);
    }

    /**
     * Asserts that the given trace correlation stats are not significantly different to
     * provided expectation value
     *
     * @param name  the name of the statistic
     * @param stats the summary of the MCMC trace of the statistic
     * @param v     the expected value for the statistic
     */
    protected void assertExpectation(String name, TraceCorrelation stats, double v) {
        double mean = stats.getMean();
        double stderr = stats.getStdErrorOfMean();
        double upper = mean + 2 * stderr;
        double lower = mean - 2 * stderr;

        assertTrue("Expected " + name + " is " + v + " but got " + mean + " +/- " + stderr,
                upper > v && lower < v);
    }
}
