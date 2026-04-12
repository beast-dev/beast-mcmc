package test.dr.inference.timeseries;

import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.core.UniformTimeGrid;
import dr.inference.timeseries.engine.gaussian.KalmanLikelihoodEngine;
import dr.inference.timeseries.gaussian.GaussianObservationModel;
import dr.inference.timeseries.gaussian.OUProcessModel;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Unit tests for {@link KalmanLikelihoodEngine}.
 * <p>
 * The engine marginalizes over the latent state and computes the log-marginal likelihood of the
 * observations via a Kalman filter forward pass.
 * <p>
 * Analytical reference values are derived for the scalar (1D state, 1D observation) case with
 * zero drift (random walk), identity design matrix, and zero-mean process.
 */
public class KalmanLikelihoodEngineTest extends TestCase {

    /** Tolerance for exact analytic comparisons. */
    private static final double TOL_EXACT = 1e-10;

    /** Tolerance for non-exact sanity checks. */
    private static final double TOL_LOOSE = 1e-6;

    public KalmanLikelihoodEngineTest(String name) {
        super(name);
    }

    // -------------------------------------------------------------------------
    // Factory helpers
    // -------------------------------------------------------------------------

    private static MatrixParameter makeMatrix(String name, double[][] values) {
        int rows = values.length;
        int cols = values[0].length;
        MatrixParameter mp = new MatrixParameter(name, rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                mp.setParameterValue(i, j, values[i][j]);
            }
        }
        return mp;
    }

    /**
     * Builds a complete scalar (1D × 1D) Kalman engine.
     *
     * @param drift     scalar drift coefficient a (A = [[a]])
     * @param diffusion scalar diffusion q (Q = [[q]])
     * @param mu0       stationary mean
     * @param p0        initial covariance (P_0 = [[p0]])
     * @param obsNoise  observation noise r (R = [[r]])
     * @param yValues   observations (one per time step; NaN = missing)
     * @param timeStep  uniform time step size
     */
    private static KalmanLikelihoodEngine makeScalar1DEngine(double drift,
                                                              double diffusion,
                                                              double mu0,
                                                              double p0,
                                                              double obsNoise,
                                                              double[] yValues,
                                                              double timeStep) {
        int T = yValues.length;
        MatrixParameter driftMat = makeMatrix("A", new double[][]{{drift}});
        MatrixParameter diffMat = makeMatrix("Q", new double[][]{{diffusion}});
        Parameter mean = new Parameter.Default(mu0);
        MatrixParameter initCov = makeMatrix("P0", new double[][]{{p0}});
        OUProcessModel process = new OUProcessModel("ou", 1, driftMat, diffMat, mean, initCov);

        MatrixParameter H = makeMatrix("H", new double[][]{{1.0}});
        MatrixParameter R = makeMatrix("R", new double[][]{{obsNoise}});
        double[][] obsData = new double[1][T];
        obsData[0] = yValues;
        MatrixParameter Y = makeMatrix("Y", obsData);
        GaussianObservationModel obs = new GaussianObservationModel("obs", 1, H, R, Y);

        TimeGrid grid = new UniformTimeGrid(T, 0.0, timeStep);
        GaussianTransitionRepresentation rep = process.getRepresentation(GaussianTransitionRepresentation.class);
        return new KalmanLikelihoodEngine(rep, obs, grid);
    }

    // -------------------------------------------------------------------------
    // Analytical test: T=1, A=0, Q=1, H=1, R=1, mu=0, P0=1, y=0
    //
    // Predicted mean = 0, predicted cov = 1
    // Innovation = 0, S = P0 + R = 2
    // LL = -0.5 * (log(2*pi) + log(2) + 0)
    // -------------------------------------------------------------------------

    public void testScalar1DExactT1ZeroObs() {
        double expected = -0.5 * (Math.log(2 * Math.PI) + Math.log(2.0));
        double actual = makeScalar1DEngine(0, 1, 0, 1, 1, new double[]{0.0}, 1.0).getLogLikelihood();
        assertEquals(expected, actual, TOL_EXACT);
    }

    // Analytical test: T=1, A=0, Q=1, H=1, R=1, mu=0, P0=1, y=5
    // Innovation = 5, S = 2
    // LL = -0.5 * (log(2*pi) + log(2) + 25/2)
    public void testScalar1DExactT1NonZeroObs() {
        double expected = -0.5 * (Math.log(2 * Math.PI) + Math.log(2.0) + 25.0 / 2.0);
        double actual = makeScalar1DEngine(0, 1, 0, 1, 1, new double[]{5.0}, 1.0).getLogLikelihood();
        assertEquals(expected, actual, TOL_EXACT);
    }

    // Analytical test: T=2, A=0, Q=1, H=1, R=1, mu=0, P0=1, y=[0,0]
    //
    // t=0: predicted=(0,1), S=2, K=0.5
    //       filtered_cov = Joseph = 0.5*1*0.5 + 0.5*1*0.5 = 0.5
    // t=1: predicted_cov = 0.5 + 1 = 1.5, S=2.5
    // LL = -0.5*(log(2*pi)+log(2)) + (-0.5*(log(2*pi)+log(2.5)))
    //    = -0.5*(2*log(2*pi) + log(5))
    public void testScalar1DExactT2ZeroObs() {
        double expected = -0.5 * (2 * Math.log(2 * Math.PI) + Math.log(5.0));
        double actual = makeScalar1DEngine(0, 1, 0, 1, 1, new double[]{0.0, 0.0}, 1.0).getLogLikelihood();
        assertEquals(expected, actual, TOL_EXACT);
    }

    // -------------------------------------------------------------------------
    // Sanity checks: finite, negative, cached
    // -------------------------------------------------------------------------

    public void testLogLikelihoodIsFinite() {
        double ll = makeScalar1DEngine(0.1, 1.0, 0.0, 1.0, 0.5,
                new double[]{1.0, -0.5, 2.0, 0.3, -1.2}, 0.1).getLogLikelihood();
        assertTrue("Log-likelihood must be finite", Double.isFinite(ll));
    }

    public void testLogLikelihoodIsNegative() {
        // For any proper probability model, the log-likelihood of a single observation is ≤ 0
        // (given a continuous density evaluated at a non-degenerate point)
        double ll = makeScalar1DEngine(0, 1, 0, 1, 1, new double[]{3.0}, 1.0).getLogLikelihood();
        assertTrue("Log-likelihood should be negative", ll < 0.0);
    }

    public void testLargerNoiseMeansHigherLikelihoodForLargeObs() {
        // With a large observation (y=10), a larger noise R gives a less negative LL
        double llSmallR = makeScalar1DEngine(0, 1, 0, 1, 1.0, new double[]{10.0}, 1.0).getLogLikelihood();
        double llLargeR = makeScalar1DEngine(0, 1, 0, 1, 100.0, new double[]{10.0}, 1.0).getLogLikelihood();
        assertTrue("Larger observation noise should increase LL for large innovations: " + llSmallR + " vs " + llLargeR,
                llLargeR > llSmallR);
    }

    // -------------------------------------------------------------------------
    // Caching
    // -------------------------------------------------------------------------

    public void testCachingReturnsSameValue() {
        KalmanLikelihoodEngine engine = makeScalar1DEngine(0, 1, 0, 1, 1, new double[]{2.0, -1.0, 0.5}, 1.0);
        double first = engine.getLogLikelihood();
        double second = engine.getLogLikelihood();
        assertEquals("Cached call must return identical value", first, second, 0.0);
    }

    public void testMakeDirtyForceRecompute() {
        KalmanLikelihoodEngine engine = makeScalar1DEngine(0, 1, 0, 1, 1, new double[]{0.0}, 1.0);
        double before = engine.getLogLikelihood();
        engine.makeDirty();
        double after = engine.getLogLikelihood();
        // Value itself should not change (no parameters changed), but recomputation must succeed
        assertEquals(before, after, TOL_EXACT);
    }

    // -------------------------------------------------------------------------
    // Missing observations
    // -------------------------------------------------------------------------

    public void testSingleMissingObservationFinite() {
        // T=1, missing → Kalman reduces to: predicted=initial, no update → LL=0
        double ll = makeScalar1DEngine(0, 1, 0, 1, 1, new double[]{Double.NaN}, 1.0).getLogLikelihood();
        assertTrue("LL with all missing obs must be finite", Double.isFinite(ll));
        assertEquals("LL with all missing obs should be 0 (no data contribution)", 0.0, ll, TOL_EXACT);
    }

    public void testMissingMiddleObservation() {
        // T=3, middle obs missing
        double ll = makeScalar1DEngine(0, 1, 0, 1, 1,
                new double[]{0.0, Double.NaN, 0.0}, 1.0).getLogLikelihood();
        assertTrue("LL with middle obs missing must be finite", Double.isFinite(ll));
        assertTrue("LL with middle obs missing must be negative", ll < 0.0);
    }

    public void testAllMissingObservations() {
        double ll = makeScalar1DEngine(0, 1, 0, 1, 1,
                new double[]{Double.NaN, Double.NaN, Double.NaN}, 1.0).getLogLikelihood();
        assertEquals("LL with all missing obs should be 0", 0.0, ll, TOL_EXACT);
    }

    // -------------------------------------------------------------------------
    // 2D state × 1D observation
    // -------------------------------------------------------------------------

    public void testTwoDimensionalStateOneDimObs() {
        int stateDim = 2;
        int obsDim = 1;
        int T = 3;

        MatrixParameter drift = makeMatrix("A", new double[][]{{0.1, 0}, {0, 0.2}});
        MatrixParameter diffusion = makeMatrix("Q", new double[][]{{1, 0}, {0, 1}});
        Parameter mean = new Parameter.Default(0.0);   // broadcast scalar
        MatrixParameter initCov = makeMatrix("P0", new double[][]{{1, 0}, {0, 1}});
        OUProcessModel process = new OUProcessModel("ou2d", stateDim, drift, diffusion, mean, initCov);

        // H = [1, 0]  (observe first state component only)
        MatrixParameter H = makeMatrix("H", new double[][]{{1.0, 0.0}});
        MatrixParameter R = makeMatrix("R", new double[][]{{1.0}});
        MatrixParameter Y = makeMatrix("Y", new double[][]{{0.5, -0.3, 1.2}});
        GaussianObservationModel obs = new GaussianObservationModel("obs", obsDim, H, R, Y);

        TimeGrid grid = new UniformTimeGrid(T, 0.0, 1.0);
        GaussianTransitionRepresentation rep = process.getRepresentation(GaussianTransitionRepresentation.class);
        KalmanLikelihoodEngine engine = new KalmanLikelihoodEngine(rep, obs, grid);

        double ll = engine.getLogLikelihood();
        assertTrue("2D state LL must be finite", Double.isFinite(ll));
        assertTrue("2D state LL must be negative", ll < 0.0);
    }

    // -------------------------------------------------------------------------
    // 2D state × 2D observation (full coverage)
    // -------------------------------------------------------------------------

    public void testTwoDimensionalFull() {
        int T = 4;

        MatrixParameter drift = makeMatrix("A", new double[][]{{0.5, 0}, {0, 0.5}});
        MatrixParameter diffusion = makeMatrix("Q", new double[][]{{1, 0}, {0, 1}});
        Parameter mean = new Parameter.Default(0.0);
        MatrixParameter initCov = makeMatrix("P0", new double[][]{{2, 0}, {0, 2}});
        OUProcessModel process = new OUProcessModel("ou2d", 2, drift, diffusion, mean, initCov);

        MatrixParameter H = makeMatrix("H", new double[][]{{1, 0}, {0, 1}});
        MatrixParameter R = makeMatrix("R", new double[][]{{0.5, 0}, {0, 0.5}});
        MatrixParameter Y = makeMatrix("Y", new double[][]{
                {1.0, -1.0, 0.5, 2.0},
                {0.0, 1.0, -0.5, -1.0}
        });
        GaussianObservationModel obs = new GaussianObservationModel("obs", 2, H, R, Y);

        TimeGrid grid = new UniformTimeGrid(T, 0.0, 0.5);
        GaussianTransitionRepresentation rep = process.getRepresentation(GaussianTransitionRepresentation.class);
        KalmanLikelihoodEngine engine = new KalmanLikelihoodEngine(rep, obs, grid);

        double ll = engine.getLogLikelihood();
        assertTrue("Full 2D LL must be finite", Double.isFinite(ll));
        assertTrue("Full 2D LL must be negative", ll < 0.0);
    }

    // -------------------------------------------------------------------------
    // Consistency: T=1 equals expected marginal normal density
    // -------------------------------------------------------------------------

    public void testScalar1DLargeInitialUncertainty() {
        // As P0 → ∞, LL should approach the flat prior limit for y at the mean
        // At least verify finite and reasonable sign
        double ll = makeScalar1DEngine(0, 1, 0, 1000.0, 1, new double[]{0.0}, 1.0).getLogLikelihood();
        assertTrue(Double.isFinite(ll));
        assertTrue(ll < 0);
    }

    // -------------------------------------------------------------------------
    // JUnit 3 boilerplate
    // -------------------------------------------------------------------------

    public static Test suite() {
        return new TestSuite(KalmanLikelihoodEngineTest.class);
    }
}
