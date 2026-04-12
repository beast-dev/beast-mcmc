package test.dr.inference.timeseries;

import dr.inference.model.MatrixParameter;
import dr.inference.timeseries.gaussian.GaussianObservationModel;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Unit tests for {@link GaussianObservationModel}.
 * <p>
 * Observations are stored as a matrix with shape (observationDimension × timeCount):
 * column j holds the observation vector at time step j.
 * A missing time step is represented by an entire column of {@link Double#NaN}.
 */
public class GaussianObservationModelTest extends TestCase {

    private static final double TOL = 1e-14;

    public GaussianObservationModelTest(String name) {
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
     * Scalar 1D model: H = [[1]], R = [[r]], y = column vector per time step.
     * observationDimension = 1, timeCount = yValues.length.
     */
    private static GaussianObservationModel makeScalar1D(double r, double... yValues) {
        MatrixParameter design = makeMatrix("H", new double[][]{{1.0}});
        MatrixParameter noise = makeMatrix("R", new double[][]{{r}});
        // observations: 1 row × timeCount cols
        double[][] obsData = new double[1][yValues.length];
        obsData[0] = yValues;
        MatrixParameter observations = makeMatrix("Y", obsData);
        return new GaussianObservationModel("obs", 1, design, noise, observations);
    }

    // -------------------------------------------------------------------------
    // Dimensions
    // -------------------------------------------------------------------------

    public void testObservationDimension() {
        GaussianObservationModel model = makeScalar1D(1.0, 0.0, 1.0, 2.0);
        assertEquals(1, model.getObservationDimension());
    }

    public void testTimeCount() {
        GaussianObservationModel model = makeScalar1D(1.0, 0.0, 1.0, 2.0);
        assertEquals(3, model.getTimeCount());
    }

    public void testTimeCount2DObs() {
        // 2D observation, 4 time steps
        MatrixParameter H = makeMatrix("H", new double[][]{{1, 0}, {0, 1}});
        MatrixParameter R = makeMatrix("R", new double[][]{{1, 0}, {0, 1}});
        MatrixParameter Y = makeMatrix("Y", new double[][]{
                {1, 2, 3, 4},   // row 0 across time
                {5, 6, 7, 8}    // row 1 across time
        });
        GaussianObservationModel model = new GaussianObservationModel("obs2d", 2, H, R, Y);
        assertEquals(2, model.getObservationDimension());
        assertEquals(4, model.getTimeCount());
    }

    // -------------------------------------------------------------------------
    // fillDesignMatrix
    // -------------------------------------------------------------------------

    public void testFillDesignMatrixIdentity1D() {
        GaussianObservationModel model = makeScalar1D(1.0, 0.0);
        double[][] H = new double[1][1];
        model.fillDesignMatrix(H);
        assertEquals(1.0, H[0][0], TOL);
    }

    public void testFillDesignMatrix2D() {
        MatrixParameter H = makeMatrix("H", new double[][]{{2.0, 3.0}, {0.0, 1.0}});
        MatrixParameter R = makeMatrix("R", new double[][]{{1, 0}, {0, 1}});
        MatrixParameter Y = makeMatrix("Y", new double[][]{{0}, {0}});
        GaussianObservationModel model = new GaussianObservationModel("obs", 2, H, R, Y);

        double[][] out = new double[2][2];
        model.fillDesignMatrix(out);
        assertEquals(2.0, out[0][0], TOL);
        assertEquals(3.0, out[0][1], TOL);
        assertEquals(0.0, out[1][0], TOL);
        assertEquals(1.0, out[1][1], TOL);
    }

    // -------------------------------------------------------------------------
    // fillNoiseCovariance
    // -------------------------------------------------------------------------

    public void testFillNoiseCovarianceScalar() {
        GaussianObservationModel model = makeScalar1D(3.5, 0.0);
        double[][] R = new double[1][1];
        model.fillNoiseCovariance(R);
        assertEquals(3.5, R[0][0], TOL);
    }

    public void testFillNoiseCovariance2D() {
        MatrixParameter H = makeMatrix("H", new double[][]{{1, 0}, {0, 1}});
        MatrixParameter R = makeMatrix("R", new double[][]{{2.0, 0.5}, {0.5, 3.0}});
        MatrixParameter Y = makeMatrix("Y", new double[][]{{0}, {0}});
        GaussianObservationModel model = new GaussianObservationModel("obs", 2, H, R, Y);

        double[][] out = new double[2][2];
        model.fillNoiseCovariance(out);
        assertEquals(2.0, out[0][0], TOL);
        assertEquals(0.5, out[0][1], TOL);
        assertEquals(0.5, out[1][0], TOL);
        assertEquals(3.0, out[1][1], TOL);
    }

    // -------------------------------------------------------------------------
    // fillObservationVector
    // -------------------------------------------------------------------------

    public void testFillObservationVector() {
        GaussianObservationModel model = makeScalar1D(1.0, 1.5, 2.5, 3.5);
        double[] obs = new double[1];

        model.fillObservationVector(0, obs);
        assertEquals(1.5, obs[0], TOL);

        model.fillObservationVector(1, obs);
        assertEquals(2.5, obs[0], TOL);

        model.fillObservationVector(2, obs);
        assertEquals(3.5, obs[0], TOL);
    }

    public void testFillObservationVector2D() {
        MatrixParameter H = makeMatrix("H", new double[][]{{1, 0}, {0, 1}});
        MatrixParameter R = makeMatrix("R", new double[][]{{1, 0}, {0, 1}});
        // t=0: (10, 20), t=1: (30, 40)
        MatrixParameter Y = makeMatrix("Y", new double[][]{
                {10.0, 30.0},
                {20.0, 40.0}
        });
        GaussianObservationModel model = new GaussianObservationModel("obs", 2, H, R, Y);

        double[] obs = new double[2];
        model.fillObservationVector(0, obs);
        assertEquals(10.0, obs[0], TOL);
        assertEquals(20.0, obs[1], TOL);

        model.fillObservationVector(1, obs);
        assertEquals(30.0, obs[0], TOL);
        assertEquals(40.0, obs[1], TOL);
    }

    public void testFillObservationVectorOutOfBoundsThrows() {
        GaussianObservationModel model = makeScalar1D(1.0, 0.0);
        try {
            model.fillObservationVector(1, new double[1]);
            fail("Expected IllegalArgumentException for out-of-bounds time index");
        } catch (IllegalArgumentException e) {
            // expected
        }
    }

    // -------------------------------------------------------------------------
    // isObservationMissing
    // -------------------------------------------------------------------------

    public void testIsObservationMissingFalse() {
        GaussianObservationModel model = makeScalar1D(1.0, 5.0);
        assertFalse(model.isObservationMissing(0));
    }

    public void testIsObservationMissingTrue() {
        GaussianObservationModel model = makeScalar1D(1.0, Double.NaN);
        assertTrue(model.isObservationMissing(0));
    }

    public void testMixedMissingAndPresent() {
        GaussianObservationModel model = makeScalar1D(1.0, 1.0, Double.NaN, 3.0);
        assertFalse(model.isObservationMissing(0));
        assertTrue(model.isObservationMissing(1));
        assertFalse(model.isObservationMissing(2));
    }

    public void testIsObservationMissing2D() {
        MatrixParameter H = makeMatrix("H", new double[][]{{1, 0}, {0, 1}});
        MatrixParameter R = makeMatrix("R", new double[][]{{1, 0}, {0, 1}});
        // t=0 full, t=1 all NaN
        MatrixParameter Y = makeMatrix("Y", new double[][]{
                {1.0, Double.NaN},
                {2.0, Double.NaN}
        });
        GaussianObservationModel model = new GaussianObservationModel("obs", 2, H, R, Y);
        assertFalse(model.isObservationMissing(0));
        assertTrue(model.isObservationMissing(1));
    }

    public void testPartialMissingThrows() {
        // Partially missing column (only one NaN in a 2D observation)
        MatrixParameter H = makeMatrix("H", new double[][]{{1, 0}, {0, 1}});
        MatrixParameter R = makeMatrix("R", new double[][]{{1, 0}, {0, 1}});
        MatrixParameter Y = makeMatrix("Y", new double[][]{
                {1.0},
                {Double.NaN}  // partial missing — should throw on access
        });
        GaussianObservationModel model = new GaussianObservationModel("obs", 2, H, R, Y);
        try {
            double[] obs = new double[2];
            model.fillObservationVector(0, obs);
            fail("Expected IllegalArgumentException for partial missing observations");
        } catch (IllegalArgumentException e) {
            // expected
        }
    }

    // -------------------------------------------------------------------------
    // Invalid construction
    // -------------------------------------------------------------------------

    public void testZeroObservationDimensionThrows() {
        try {
            MatrixParameter H = makeMatrix("H", new double[][]{{1}});
            MatrixParameter R = makeMatrix("R", new double[][]{{1}});
            MatrixParameter Y = makeMatrix("Y", new double[][]{{0}});
            new GaussianObservationModel("bad", 0, H, R, Y);
            fail("Expected IllegalArgumentException for observationDimension = 0");
        } catch (IllegalArgumentException e) {
            // expected
        }
    }

    public void testNullDesignMatrixThrows() {
        try {
            MatrixParameter R = makeMatrix("R", new double[][]{{1}});
            MatrixParameter Y = makeMatrix("Y", new double[][]{{0}});
            new GaussianObservationModel("bad", 1, null, R, Y);
            fail("Expected IllegalArgumentException for null design matrix");
        } catch (IllegalArgumentException e) {
            // expected
        }
    }

    // -------------------------------------------------------------------------
    // JUnit 3 boilerplate
    // -------------------------------------------------------------------------

    public static Test suite() {
        return new TestSuite(GaussianObservationModelTest.class);
    }
}
