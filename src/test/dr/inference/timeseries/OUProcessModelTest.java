package test.dr.inference.timeseries;

import static test.dr.inference.timeseries.OUTimeSeriesTestSupport.representation;
import static test.dr.inference.timeseries.OUTimeSeriesTestSupport.representable;
import static test.dr.inference.timeseries.OUTimeSeriesTestSupport.latent;
import static test.dr.inference.timeseries.OUTimeSeriesTestSupport.supportsRepresentation;
import static test.dr.inference.timeseries.OUTimeSeriesTestSupport.getTransitionMatrix;
import static test.dr.inference.timeseries.OUTimeSeriesTestSupport.getTransitionOffset;
import static test.dr.inference.timeseries.OUTimeSeriesTestSupport.getTransitionCovariance;

import dr.inference.model.MatrixParameter;
import dr.inference.model.GivensRotationMatrixParameter;
import dr.inference.model.OrthogonalBlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.core.UniformTimeGrid;
import dr.inference.timeseries.gaussian.EulerOUProcessModel;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Unit tests for {@link OUProcessModel} (exact matrix-exponential transitions) and
 * {@link EulerOUProcessModel} (first-order Euler approximation).
 *
 * <p>Exact reference values for the scalar OU process are:
 * <pre>
 *   F(dt) = exp(-a dt)
 *   b(dt) = (1 - exp(-a dt)) mu
 *   V(dt) = (q / 2a) * (1 - exp(-2a dt))    [1-D, a != 0]
 *         = q * dt                            [a = 0]
 * </pre>
 */
public class OUProcessModelTest extends TestCase {

    /** Tight tolerance for quantities that do not pass through expm. */
    private static final double TOL = 1e-14;

    /**
     * Tolerance for values produced by the Padé-based matrix exponential; the
     * scaling-and-squaring algorithm is accurate to near machine precision but a
     * slightly relaxed threshold avoids brittle last-bit comparisons.
     */
    private static final double EXACT_TOL = 1e-12;

    public OUProcessModelTest(String name) {
        super(name);
    }

    // ── Factory helpers ──────────────────────────────────────────────────────────

    private static MatrixParameter makeMatrix(String name, double[][] values) {
        MatrixParameter mp = new MatrixParameter(name, values.length, values[0].length);
        for (int i = 0; i < values.length; i++) {
            for (int j = 0; j < values[0].length; j++) {
                mp.setParameterValue(i, j, values[i][j]);
            }
        }
        return mp;
    }

    private static OUProcessModel makeScalar1D(double a, double q, double mu0, double p0) {
        return new OUProcessModel("ou1d", 1,
                makeMatrix("A", new double[][]{{a}}),
                makeMatrix("Q", new double[][]{{q}}),
                new Parameter.Default(mu0),
                makeMatrix("P0", new double[][]{{p0}}));
    }

    private static OUProcessModel make2D(double a11, double a22,
                                         double q11, double q22,
                                         double mu0, double p0Diag) {
        return new OUProcessModel("ou2d", 2,
                makeMatrix("A", new double[][]{{a11, 0}, {0, a22}}),
                makeMatrix("Q", new double[][]{{q11, 0}, {0, q22}}),
                new Parameter.Default(mu0),
                makeMatrix("P0", new double[][]{{p0Diag, 0}, {0, p0Diag}}));
    }

    private static OrthogonalBlockDiagonalPolarStableMatrixParameter makeOrthogonalBlock2D(final String name) {
        final Parameter angles = new Parameter.Default(name + ".angles", new double[]{0.2});
        final GivensRotationMatrixParameter rotation =
                new GivensRotationMatrixParameter(name + ".rotation", 2, angles);
        final Parameter scalar = new Parameter.Default(0);
        final Parameter rho = new Parameter.Default(name + ".rho", new double[]{0.85});
        final Parameter theta = new Parameter.Default(name + ".theta", new double[]{0.25});
        final Parameter t = new Parameter.Default(name + ".t", new double[]{-0.08});
        return new OrthogonalBlockDiagonalPolarStableMatrixParameter(
                name, rotation, scalar, rho, theta, t);
    }

    private static EulerOUProcessModel makeEulerScalar1D(double a, double q, double mu0, double p0) {
        return new EulerOUProcessModel("euler1d", 1,
                makeMatrix("A", new double[][]{{a}}),
                makeMatrix("Q", new double[][]{{q}}),
                new Parameter.Default(mu0),
                makeMatrix("P0", new double[][]{{p0}}));
    }

    private static EulerOUProcessModel makeEuler2D(double a11, double a22,
                                                    double q11, double q22,
                                                    double mu0, double p0Diag) {
        return new EulerOUProcessModel("euler2d", 2,
                makeMatrix("A", new double[][]{{a11, 0}, {0, a22}}),
                makeMatrix("Q", new double[][]{{q11, 0}, {0, q22}}),
                new Parameter.Default(mu0),
                makeMatrix("P0", new double[][]{{p0Diag, 0}, {0, p0Diag}}));
    }

    // ── Dimension ────────────────────────────────────────────────────────────────

    public void testGetStateDimension1D() {
        assertEquals(1, makeScalar1D(0.5, 1.0, 0.0, 1.0).getStateDimension());
    }

    public void testGetStateDimension2D() {
        assertEquals(2, make2D(1, 1, 1, 1, 0, 1).getStateDimension());
    }

    public void testDefaultCovarianceGradientMethodUsesVanLoanForDenseDrift() {
        final OUProcessModel process = make2D(1.0, 2.0, 1.5, 1.7, 0.0, 1.0);
        assertEquals(OUProcessModel.CovarianceGradientMethod.VAN_LOAN_ADJOINT,
                process.getCovarianceGradientMethod());
    }

    public void testDefaultCovarianceGradientMethodUsesStationaryLyapunovForOrthogonalBlockDrift() {
        final OUProcessModel process = new OUProcessModel(
                "ou2dBlockDefault",
                2,
                makeOrthogonalBlock2D("Ablock.default"),
                makeMatrix("Q.default", new double[][]{{1.3, 0.1}, {0.1, 1.7}}),
                new Parameter.Default(0.0),
                makeMatrix("P0.default", new double[][]{{1.0, 0.0}, {0.0, 1.0}}));
        assertEquals(OUProcessModel.CovarianceGradientMethod.STATIONARY_LYAPUNOV,
                process.getCovarianceGradientMethod());
    }

    public void testUsesOrthogonalBlockSelectionChartDetectsExpectedCharts() {
        assertFalse(OUProcessModel.usesOrthogonalBlockSelectionChart(
                makeMatrix("A.dense", new double[][]{{1.0, -0.2}, {0.1, 1.4}})));
        assertTrue(OUProcessModel.usesOrthogonalBlockSelectionChart(
                makeOrthogonalBlock2D("A.block.chart")));
    }

    // ── Initial mean ─────────────────────────────────────────────────────────────

    public void testInitialMeanScalarBroadcast() {
        double[] mean = new double[2];
        make2D(0, 0, 1, 1, 3.0, 1.0).getInitialMean(mean);
        assertEquals(3.0, mean[0], TOL);
        assertEquals(3.0, mean[1], TOL);
    }

    public void testInitialMeanVector() {
        OUProcessModel ou = new OUProcessModel("ouVec", 2,
                makeMatrix("A", new double[][]{{0, 0}, {0, 0}}),
                makeMatrix("Q", new double[][]{{1, 0}, {0, 1}}),
                new Parameter.Default(new double[]{1.5, 2.5}),
                makeMatrix("P0", new double[][]{{1, 0}, {0, 1}}));
        double[] mean = new double[2];
        ou.getInitialMean(mean);
        assertEquals(1.5, mean[0], TOL);
        assertEquals(2.5, mean[1], TOL);
    }

    public void testInitialMean1D() {
        double[] mean = new double[1];
        makeScalar1D(0.0, 1.0, 7.0, 1.0).getInitialMean(mean);
        assertEquals(7.0, mean[0], TOL);
    }

    // ── Initial covariance ───────────────────────────────────────────────────────

    public void testInitialCovariance1D() {
        double[][] cov = new double[1][1];
        makeScalar1D(0.0, 1.0, 0.0, 4.0).getInitialCovariance(cov);
        assertEquals(4.0, cov[0][0], TOL);
    }

    public void testInitialCovariance2D() {
        double[][] cov = new double[2][2];
        make2D(0, 0, 1, 1, 0, 2.0).getInitialCovariance(cov);
        assertEquals(2.0, cov[0][0], TOL);
        assertEquals(2.0, cov[1][1], TOL);
        assertEquals(0.0, cov[0][1], TOL);
        assertEquals(0.0, cov[1][0], TOL);
    }

    // ── Exact transition matrix: F = expm(-A dt) ─────────────────────────────────

    public void testTransitionMatrixZeroDrift() {
        // A = 0 → expm(0) = I
        TimeGrid grid = new UniformTimeGrid(2, 0.0, 1.0);
        double[][] F = new double[1][1];
        getTransitionMatrix(makeScalar1D(0.0, 1.0, 0.0, 1.0), 0, 1, grid, F);
        assertEquals(1.0, F[0][0], TOL);
    }

    public void testTransitionMatrix1DWithDrift() {
        // a = 2, dt = 0.5  →  F = exp(-1)
        TimeGrid grid = new UniformTimeGrid(2, 0.0, 0.5);
        double[][] F = new double[1][1];
        getTransitionMatrix(makeScalar1D(2.0, 1.0, 0.0, 1.0), 0, 1, grid, F);
        assertEquals(Math.exp(-1.0), F[0][0], EXACT_TOL);
    }

    public void testTransitionMatrix2DDiagonal() {
        // A = diag(1, 2), dt = 1  →  F = diag(exp(-1), exp(-2))
        TimeGrid grid = new UniformTimeGrid(2, 0.0, 1.0);
        double[][] F = new double[2][2];
        getTransitionMatrix(make2D(1.0, 2.0, 1.0, 1.0, 0.0, 1.0), 0, 1, grid, F);
        assertEquals(Math.exp(-1.0), F[0][0], EXACT_TOL);
        assertEquals(Math.exp(-2.0), F[1][1], EXACT_TOL);
        assertEquals(0.0, F[0][1], EXACT_TOL);
        assertEquals(0.0, F[1][0], EXACT_TOL);
    }

    public void testTransitionMatrixNonUnitStep() {
        // a = 1, dt = 0.25  →  F = exp(-0.25)
        TimeGrid grid = new UniformTimeGrid(3, 0.0, 0.25);
        double[][] F = new double[1][1];
        getTransitionMatrix(makeScalar1D(1.0, 1.0, 0.0, 1.0), 1, 2, grid, F);
        assertEquals(Math.exp(-0.25), F[0][0], EXACT_TOL);
    }

    // ── Exact transition offset: b = (I - F) mu ──────────────────────────────────

    public void testTransitionOffsetZeroMean() {
        // mu = 0 → b = 0 regardless of A
        TimeGrid grid = new UniformTimeGrid(2, 0.0, 1.0);
        double[] offset = new double[1];
        getTransitionOffset(makeScalar1D(2.0, 1.0, 0.0, 1.0), 0, 1, grid, offset);
        assertEquals(0.0, offset[0], EXACT_TOL);
    }

    public void testTransitionOffsetZeroDrift() {
        // A = 0 → F = I → b = (I - I) mu = 0
        TimeGrid grid = new UniformTimeGrid(2, 0.0, 1.0);
        double[] offset = new double[1];
        getTransitionOffset(makeScalar1D(0.0, 1.0, 5.0, 1.0), 0, 1, grid, offset);
        assertEquals(0.0, offset[0], EXACT_TOL);
    }

    public void testTransitionOffset1D() {
        // a = 2, mu = 3, dt = 0.5  →  b = (1 - exp(-1)) * 3
        TimeGrid grid = new UniformTimeGrid(2, 0.0, 0.5);
        double[] offset = new double[1];
        getTransitionOffset(makeScalar1D(2.0, 1.0, 3.0, 1.0), 0, 1, grid, offset);
        assertEquals(3.0 * (1.0 - Math.exp(-1.0)), offset[0], EXACT_TOL);
    }

    // ── Exact transition covariance (Van Loan) ───────────────────────────────────

    public void testTransitionCovariance1DZeroDrift() {
        // A = 0 → V = Q * dt  (exact integral reduces to Euler)
        TimeGrid grid = new UniformTimeGrid(2, 0.0, 2.0);
        double[][] Qt = new double[1][1];
        getTransitionCovariance(makeScalar1D(0.0, 3.0, 0.0, 1.0), 0, 1, grid, Qt);
        assertEquals(6.0, Qt[0][0], EXACT_TOL);
    }

    public void testTransitionCovariance1DWithDrift() {
        // a = 1, q = 2, dt = 1  →  V = (q/2a)(1 - exp(-2a dt)) = 1 * (1 - exp(-2))
        TimeGrid grid = new UniformTimeGrid(2, 0.0, 1.0);
        double[][] Qt = new double[1][1];
        getTransitionCovariance(makeScalar1D(1.0, 2.0, 0.0, 1.0), 0, 1, grid, Qt);
        final double expected = (2.0 / 2.0) * (1.0 - Math.exp(-2.0));
        assertEquals(expected, Qt[0][0], EXACT_TOL);
    }

    public void testTransitionCovariance2DZeroDrift() {
        // A = 0 → V = Q * dt
        TimeGrid grid = new UniformTimeGrid(2, 0.0, 0.5);
        double[][] Qt = new double[2][2];
        getTransitionCovariance(make2D(0, 0, 2.0, 4.0, 0, 1.0), 0, 1, grid, Qt);
        assertEquals(1.0, Qt[0][0], EXACT_TOL);
        assertEquals(2.0, Qt[1][1], EXACT_TOL);
        assertEquals(0.0, Qt[0][1], EXACT_TOL);
        assertEquals(0.0, Qt[1][0], EXACT_TOL);
    }

    // ── Time-series adapter bridge ────────────────────────────────────────────────

    public void testSupportsGaussianRepresentation() {
        assertTrue(supportsRepresentation(
                makeScalar1D(0, 1, 0, 1),
                GaussianTransitionRepresentation.class));
    }

    public void testGetRepresentationReturnsTransitionRepresentation() {
        OUProcessModel ou = makeScalar1D(0, 1, 0, 1);
        GaussianTransitionRepresentation representation =
                representation(ou, GaussianTransitionRepresentation.class);
        assertNotNull(representation);
        assertFalse(representation == ou);
    }

    public void testGetRepresentationUnsupportedThrows() {
        try {
            representation(makeScalar1D(0, 1, 0, 1), String.class);
            fail("Expected IllegalArgumentException");
        } catch (IllegalArgumentException e) { /* expected */ }
    }

    // ── Invalid construction ─────────────────────────────────────────────────────

    public void testZeroStateDimensionThrows() {
        try {
            new OUProcessModel("bad", 0,
                    makeMatrix("A", new double[][]{{1}}),
                    makeMatrix("Q", new double[][]{{1}}),
                    new Parameter.Default(0.0),
                    makeMatrix("P0", new double[][]{{1}}));
            fail("Expected IllegalArgumentException");
        } catch (IllegalArgumentException e) { /* expected */ }
    }

    public void testWrongDriftShapeThrows() {
        try {
            new OUProcessModel("bad", 2,
                    makeMatrix("A", new double[][]{{1}}),
                    makeMatrix("Q", new double[][]{{1, 0}, {0, 1}}),
                    new Parameter.Default(0.0),
                    makeMatrix("P0", new double[][]{{1, 0}, {0, 1}}));
            fail("Expected IllegalArgumentException");
        } catch (IllegalArgumentException e) { /* expected */ }
    }

    public void testWrongMeanDimensionThrows() {
        try {
            new OUProcessModel("bad", 2,
                    makeMatrix("A", new double[][]{{1, 0}, {0, 1}}),
                    makeMatrix("Q", new double[][]{{1, 0}, {0, 1}}),
                    new Parameter.Default(new double[]{1.0, 2.0, 3.0}),
                    makeMatrix("P0", new double[][]{{1, 0}, {0, 1}}));
            fail("Expected IllegalArgumentException");
        } catch (IllegalArgumentException e) { /* expected */ }
    }

    // ── EulerOUProcessModel: verify Euler-specific values ────────────────────────
    // These tests document and guard the first-order approximations.

    public void testEulerTransitionMatrix1DWithDrift() {
        // F = I - dt * A = 1 - 0.5 * 2 = 0
        TimeGrid grid = new UniformTimeGrid(2, 0.0, 0.5);
        double[][] F = new double[1][1];
        getTransitionMatrix(makeEulerScalar1D(2.0, 1.0, 0.0, 1.0), 0, 1, grid, F);
        assertEquals(0.0, F[0][0], TOL);
    }

    public void testEulerTransitionMatrix2DDiagonal() {
        // F = I - dt * diag(1, 2) = diag(0, -1)
        TimeGrid grid = new UniformTimeGrid(2, 0.0, 1.0);
        double[][] F = new double[2][2];
        getTransitionMatrix(makeEuler2D(1.0, 2.0, 1.0, 1.0, 0.0, 1.0), 0, 1, grid, F);
        assertEquals(0.0, F[0][0], TOL);
        assertEquals(-1.0, F[1][1], TOL);
        assertEquals(0.0, F[0][1], TOL);
        assertEquals(0.0, F[1][0], TOL);
    }

    public void testEulerTransitionOffset1D() {
        // f = dt * A * mu = 0.5 * 2 * 3 = 3
        TimeGrid grid = new UniformTimeGrid(2, 0.0, 0.5);
        double[] offset = new double[1];
        getTransitionOffset(makeEulerScalar1D(2.0, 1.0, 3.0, 1.0), 0, 1, grid, offset);
        assertEquals(3.0, offset[0], TOL);
    }

    public void testEulerTransitionCovariance1D() {
        // V = dt * Q = 2 * 3 = 6
        TimeGrid grid = new UniformTimeGrid(2, 0.0, 2.0);
        double[][] Qt = new double[1][1];
        getTransitionCovariance(makeEulerScalar1D(0.0, 3.0, 0.0, 1.0), 0, 1, grid, Qt);
        assertEquals(6.0, Qt[0][0], TOL);
    }

    // ── JUnit 3 boilerplate ──────────────────────────────────────────────────────

    public static Test suite() {
        return new TestSuite(OUProcessModelTest.class);
    }
}
