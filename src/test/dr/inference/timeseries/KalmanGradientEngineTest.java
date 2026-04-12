package test.dr.inference.timeseries;

import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.core.UniformTimeGrid;
import dr.inference.timeseries.engine.gaussian.KalmanGradientEngine;
import dr.inference.timeseries.engine.gaussian.KalmanLikelihoodEngine;
import dr.inference.timeseries.gaussian.GaussianObservationModel;
import dr.inference.timeseries.gaussian.OUProcessModel;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Unit tests for {@link KalmanGradientEngine}.
 * <p>
 * The engine uses central-difference numerical differentiation.  Tests verify:
 * <ul>
 *   <li>Which parameters are supported.</li>
 *   <li>Gradient elements are finite.</li>
 *   <li>Gradient values agree with an independent external finite-difference check.</li>
 *   <li>Analytic gradient values for the simplest 1D scalar case.</li>
 * </ul>
 */
public class KalmanGradientEngineTest extends TestCase {

    /** Tolerance for comparing with externally-computed finite-difference gradients. */
    private static final double TOL_FD = 1e-5;

    /** Tolerance for analytic gradient comparisons. */
    private static final double TOL_ANALYTIC = 1e-6;

    public KalmanGradientEngineTest(String name) {
        super(name);
    }

    // -------------------------------------------------------------------------
    // Inner container for the assembled model components
    // -------------------------------------------------------------------------

    private static class Model1D {
        final OUProcessModel process;
        final GaussianObservationModel obs;
        final KalmanLikelihoodEngine likelihoodEngine;
        final KalmanGradientEngine gradientEngine;

        Model1D(OUProcessModel process,
                GaussianObservationModel obs,
                KalmanLikelihoodEngine likelihoodEngine,
                KalmanGradientEngine gradientEngine) {
            this.process = process;
            this.obs = obs;
            this.likelihoodEngine = likelihoodEngine;
            this.gradientEngine = gradientEngine;
        }
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
     * Builds a scalar 1D engine for gradient testing with a non-zero observation so all gradients
     * are non-degenerate (y=5 gives a non-trivial quadratic term in the LL).
     */
    private static Model1D makeScalar1D(double a, double q, double mu0, double p0,
                                         double r, double[] yValues, double dt) {
        int T = yValues.length;

        MatrixParameter drift = makeMatrix("A", new double[][]{{a}});
        MatrixParameter diffusion = makeMatrix("Q", new double[][]{{q}});
        Parameter mean = new Parameter.Default(mu0);
        MatrixParameter initCov = makeMatrix("P0", new double[][]{{p0}});
        OUProcessModel process = new OUProcessModel("ou", 1, drift, diffusion, mean, initCov);

        MatrixParameter H = makeMatrix("H", new double[][]{{1.0}});
        MatrixParameter R = makeMatrix("R", new double[][]{{r}});
        double[][] obsData = new double[1][T];
        obsData[0] = yValues;
        MatrixParameter Y = makeMatrix("Y", obsData);
        GaussianObservationModel obs = new GaussianObservationModel("obs", 1, H, R, Y);

        TimeGrid grid = new UniformTimeGrid(T, 0.0, dt);
        GaussianTransitionRepresentation rep = process.getRepresentation(GaussianTransitionRepresentation.class);
        KalmanLikelihoodEngine likelihoodEngine = new KalmanLikelihoodEngine(rep, obs, grid);
        KalmanGradientEngine gradientEngine = new KalmanGradientEngine(likelihoodEngine, process, obs);

        return new Model1D(process, obs, likelihoodEngine, gradientEngine);
    }

    // -------------------------------------------------------------------------
    // supportsGradientWrt
    // -------------------------------------------------------------------------

    public void testSupportsGradientWrtDriftMatrix() {
        Model1D m = makeScalar1D(0, 1, 0, 1, 1, new double[]{5.0}, 1.0);
        assertTrue(m.gradientEngine.supportsGradientWrt(m.process.getDriftMatrix()));
    }

    public void testSupportsGradientWrtDiffusionMatrix() {
        Model1D m = makeScalar1D(0, 1, 0, 1, 1, new double[]{5.0}, 1.0);
        assertTrue(m.gradientEngine.supportsGradientWrt(m.process.getDiffusionMatrix()));
    }

    public void testSupportsGradientWrtStationaryMean() {
        Model1D m = makeScalar1D(0, 1, 0, 1, 1, new double[]{5.0}, 1.0);
        assertTrue(m.gradientEngine.supportsGradientWrt(m.process.getStationaryMeanParameter()));
    }

    public void testSupportsGradientWrtInitialCovariance() {
        Model1D m = makeScalar1D(0, 1, 0, 1, 1, new double[]{5.0}, 1.0);
        assertTrue(m.gradientEngine.supportsGradientWrt(m.process.getInitialCovarianceParameter()));
    }

    public void testSupportsGradientWrtDesignMatrix() {
        Model1D m = makeScalar1D(0, 1, 0, 1, 1, new double[]{5.0}, 1.0);
        assertTrue(m.gradientEngine.supportsGradientWrt(m.obs.getDesignMatrix()));
    }

    public void testSupportsGradientWrtNoiseCovariance() {
        Model1D m = makeScalar1D(0, 1, 0, 1, 1, new double[]{5.0}, 1.0);
        assertTrue(m.gradientEngine.supportsGradientWrt(m.obs.getNoiseCovariance()));
    }

    public void testSupportsGradientWrtObservations() {
        Model1D m = makeScalar1D(0, 1, 0, 1, 1, new double[]{5.0}, 1.0);
        assertTrue(m.gradientEngine.supportsGradientWrt(m.obs.getObservations()));
    }

    public void testSupportsGradientWrtNull() {
        Model1D m = makeScalar1D(0, 1, 0, 1, 1, new double[]{5.0}, 1.0);
        assertFalse(m.gradientEngine.supportsGradientWrt(null));
    }

    public void testSupportsGradientWrtUnknownParameter() {
        Model1D m = makeScalar1D(0, 1, 0, 1, 1, new double[]{5.0}, 1.0);
        Parameter unrelated = new Parameter.Default(42.0);
        assertFalse(m.gradientEngine.supportsGradientWrt(unrelated));
    }

    public void testUnsupportedParameterThrows() {
        Model1D m = makeScalar1D(0, 1, 0, 1, 1, new double[]{5.0}, 1.0);
        try {
            m.gradientEngine.getGradientWrt(new Parameter.Default(1.0));
            fail("Expected IllegalArgumentException for unsupported parameter");
        } catch (IllegalArgumentException e) {
            // expected
        }
    }

    // -------------------------------------------------------------------------
    // Gradient is finite
    // -------------------------------------------------------------------------

    public void testGradientWrtNoiseCovarianceIsFinite() {
        Model1D m = makeScalar1D(0, 1, 0, 1, 1, new double[]{5.0}, 1.0);
        double[] grad = m.gradientEngine.getGradientWrt(m.obs.getNoiseCovariance());
        assertEquals(1, grad.length);
        assertTrue("Gradient wrt R must be finite", Double.isFinite(grad[0]));
    }

    public void testGradientWrtDriftMatrixIsFinite() {
        Model1D m = makeScalar1D(0.1, 1, 0, 1, 1, new double[]{1.0, -0.5}, 0.5);
        double[] grad = m.gradientEngine.getGradientWrt(m.process.getDriftMatrix());
        assertEquals(1, grad.length);
        assertTrue("Gradient wrt A must be finite", Double.isFinite(grad[0]));
    }

    public void testGradientWrtAllParametersFinite() {
        Model1D m = makeScalar1D(0.2, 1.5, 0.5, 2.0, 0.8, new double[]{1.0, -1.0, 2.0}, 0.25);

        Parameter[] params = {
                m.process.getDriftMatrix(),
                m.process.getDiffusionMatrix(),
                m.process.getStationaryMeanParameter(),
                m.process.getInitialCovarianceParameter(),
                m.obs.getDesignMatrix(),
                m.obs.getNoiseCovariance(),
                m.obs.getObservations()
        };

        for (Parameter p : params) {
            double[] grad = m.gradientEngine.getGradientWrt(p);
            for (int i = 0; i < grad.length; i++) {
                assertTrue("Gradient element " + i + " of " + p.getId() + " must be finite",
                        Double.isFinite(grad[i]));
            }
        }
    }

    // -------------------------------------------------------------------------
    // Analytic gradient check for R (scalar 1D, T=1)
    //
    // LL = -0.5 * (log(2*pi) + log(S) + e^2/S),  S = P0 + R,  e = y - mu
    // d(LL)/dR = -0.5 * (1/S - e^2/S^2)
    //
    // With y=5, mu=0, P0=1, R=1: S=2, e=5
    // d(LL)/dR = -0.5*(0.5 - 25/4) = -0.5*(-5.75) = 2.875
    // -------------------------------------------------------------------------

    public void testAnalyticGradientWrtR_T1() {
        Model1D m = makeScalar1D(0, 1, 0, 1, 1.0, new double[]{5.0}, 1.0);
        double[] grad = m.gradientEngine.getGradientWrt(m.obs.getNoiseCovariance());
        assertEquals(1, grad.length);
        assertEquals("d(LL)/dR should match analytic value", 2.875, grad[0], TOL_ANALYTIC);
    }

    // Gradient wrt y (observation): d(LL)/dy = -e/S = -(y-mu)/S = -5/2 = -2.5
    public void testAnalyticGradientWrtY_T1() {
        Model1D m = makeScalar1D(0, 1, 0, 1, 1.0, new double[]{5.0}, 1.0);
        double[] grad = m.gradientEngine.getGradientWrt(m.obs.getObservations());
        assertEquals(1, grad.length);
        assertEquals("d(LL)/dy should match analytic value", -2.5, grad[0], TOL_ANALYTIC);
    }

    // Gradient wrt P0 (initial covariance): d(LL)/dP0 = -0.5*(1/S - e^2/S^2)
    // same as d/dR since S = P0 + R
    public void testAnalyticGradientWrtP0_T1() {
        Model1D m = makeScalar1D(0, 1, 0, 1.0, 1.0, new double[]{5.0}, 1.0);
        double[] grad = m.gradientEngine.getGradientWrt(m.process.getInitialCovarianceParameter());
        assertEquals(1, grad.length);
        assertEquals("d(LL)/dP0 should match analytic value", 2.875, grad[0], TOL_ANALYTIC);
    }

    // -------------------------------------------------------------------------
    // Independent finite-difference verification
    //
    // Compute gradient externally by perturbing each parameter and comparing
    // with what the engine returns.
    // -------------------------------------------------------------------------

    public void testGradientMatchesExternalFiniteDifference() {
        double a = 0.1, q = 1.5, mu = 0.5, p0 = 2.0, r = 0.8;
        double[] yValues = {1.0, -0.5, 2.0};
        double dt = 0.5;

        Model1D m = makeScalar1D(a, q, mu, p0, r, yValues, dt);

        // Check gradient wrt noiseCovariance via an external central-difference
        double h = 1e-6;
        Parameter R = m.obs.getNoiseCovariance();

        double origR = R.getParameterValue(0);
        R.setParameterValue(0, origR + h);
        m.likelihoodEngine.makeDirty();
        double llPlus = m.likelihoodEngine.getLogLikelihood();

        R.setParameterValue(0, origR - h);
        m.likelihoodEngine.makeDirty();
        double llMinus = m.likelihoodEngine.getLogLikelihood();

        R.setParameterValue(0, origR);
        m.likelihoodEngine.makeDirty();

        double externalFD = (llPlus - llMinus) / (2 * h);
        double[] engineGrad = m.gradientEngine.getGradientWrt(R);

        assertEquals("Engine gradient wrt R must match external FD",
                externalFD, engineGrad[0], TOL_FD);
    }

    public void testGradientWrtDriftMatchesExternalFD() {
        double a = 0.3, q = 1.0, mu = 1.0, p0 = 1.0, r = 1.0;
        double[] yValues = {2.0, 1.5, -0.5, 3.0};
        double dt = 0.25;

        Model1D m = makeScalar1D(a, q, mu, p0, r, yValues, dt);

        double h = 1e-6;
        Parameter A = m.process.getDriftMatrix();

        double origA = A.getParameterValue(0);
        A.setParameterValue(0, origA + h);
        m.likelihoodEngine.makeDirty();
        double llPlus = m.likelihoodEngine.getLogLikelihood();

        A.setParameterValue(0, origA - h);
        m.likelihoodEngine.makeDirty();
        double llMinus = m.likelihoodEngine.getLogLikelihood();

        A.setParameterValue(0, origA);
        m.likelihoodEngine.makeDirty();

        double externalFD = (llPlus - llMinus) / (2 * h);
        double[] engineGrad = m.gradientEngine.getGradientWrt(A);

        assertEquals("Engine gradient wrt A must match external FD",
                externalFD, engineGrad[0], TOL_FD);
    }

    // -------------------------------------------------------------------------
    // JUnit 3 boilerplate
    // -------------------------------------------------------------------------

    public static Test suite() {
        return new TestSuite(KalmanGradientEngineTest.class);
    }
}
