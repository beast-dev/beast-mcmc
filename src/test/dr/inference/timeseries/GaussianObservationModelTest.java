package test.dr.inference.timeseries;

import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.timeseries.model.gaussian.LinearGaussianObservationModel;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Unit tests for {@link LinearGaussianObservationModel}.
 * <p>
 * Observations are stored as a matrix with shape (observationDimension × timeCount):
 * column j holds the observation vector at time step j.
 * A fully missing time step is an entire column of {@link Double#NaN}.
 * A partially observed column is accepted by the model (the canonical engines mask it per
 * cell) but rejected by the complete-column accessor used by the moment-form debug engines.
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

    private static MatrixParameter identity2() {
        return makeMatrix("I2", new double[][]{{1, 0}, {0, 1}});
    }

    /**
     * Scalar 1D model: H = [[1]], R = [[r]], y = column vector per time step.
     * observationDimension = 1, timeCount = yValues.length.
     */
    private static LinearGaussianObservationModel makeScalar1D(double r, double... yValues) {
        MatrixParameter design = makeMatrix("H", new double[][]{{1.0}});
        MatrixParameter noise = makeMatrix("R", new double[][]{{r}});
        // observations: 1 row × timeCount cols
        double[][] obsData = new double[1][yValues.length];
        obsData[0] = yValues;
        MatrixParameter observations = makeMatrix("Y", obsData);
        return new LinearGaussianObservationModel("obs", 1, design, noise, observations);
    }

    private static IllegalArgumentException expectConstructionFailure(int observationDimension,
                                                                      MatrixParameter design,
                                                                      MatrixParameter noise,
                                                                      MatrixParameter observations) {
        try {
            new LinearGaussianObservationModel("bad", observationDimension, design, noise, observations);
        } catch (IllegalArgumentException e) {
            return e;
        }
        fail("Expected IllegalArgumentException from the constructor");
        return null;
    }

    // -------------------------------------------------------------------------
    // Dimensions
    // -------------------------------------------------------------------------

    public void testObservationDimension() {
        LinearGaussianObservationModel model = makeScalar1D(1.0, 0.0, 1.0, 2.0);
        assertEquals(1, model.getObservationDimension());
    }

    public void testTimeCount() {
        LinearGaussianObservationModel model = makeScalar1D(1.0, 0.0, 1.0, 2.0);
        assertEquals(3, model.getTimeCount());
    }

    public void testTimeCount2DObs() {
        // 2D observation, 4 time steps
        MatrixParameter H = identity2();
        MatrixParameter R = identity2();
        MatrixParameter Y = makeMatrix("Y", new double[][]{
                {1, 2, 3, 4},   // row 0 across time
                {5, 6, 7, 8}    // row 1 across time
        });
        LinearGaussianObservationModel model = new LinearGaussianObservationModel("obs2d", 2, H, R, Y);
        assertEquals(2, model.getObservationDimension());
        assertEquals(4, model.getTimeCount());
    }

    // -------------------------------------------------------------------------
    // fillDesignMatrix
    // -------------------------------------------------------------------------

    public void testFillDesignMatrixIdentity1D() {
        LinearGaussianObservationModel model = makeScalar1D(1.0, 0.0);
        double[] H = new double[1];
        model.fillDesignMatrixFlat(H, 1);
        assertEquals(1.0, H[0], TOL);
    }

    public void testFillDesignMatrix2D() {
        MatrixParameter H = makeMatrix("H", new double[][]{{2.0, 3.0}, {0.0, 1.0}});
        MatrixParameter R = identity2();
        MatrixParameter Y = makeMatrix("Y", new double[][]{{0}, {0}});
        LinearGaussianObservationModel model = new LinearGaussianObservationModel("obs", 2, H, R, Y);

        double[] out = new double[4];
        model.fillDesignMatrixFlat(out, 2);
        assertEquals(2.0, out[0], TOL);
        assertEquals(3.0, out[1], TOL);
        assertEquals(0.0, out[2], TOL);
        assertEquals(1.0, out[3], TOL);
    }

    public void testFillRectangularDesignMatrix() {
        // 2 observed coordinates of a 3-dimensional state.
        MatrixParameter H = makeMatrix("H", new double[][]{{1.0, 0.0, 0.5}, {0.0, 1.0, -0.5}});
        MatrixParameter R = identity2();
        MatrixParameter Y = makeMatrix("Y", new double[][]{{0}, {0}});
        LinearGaussianObservationModel model = new LinearGaussianObservationModel("obs", 2, H, R, Y);

        double[] out = new double[6];
        model.fillDesignMatrixFlat(out, 3);
        assertEquals(0.5, out[2], TOL);
        assertEquals(-0.5, out[5], TOL);
    }

    // -------------------------------------------------------------------------
    // fillNoiseCovariance
    // -------------------------------------------------------------------------

    public void testFillNoiseCovarianceScalar() {
        LinearGaussianObservationModel model = makeScalar1D(3.5, 0.0);
        double[] R = new double[1];
        model.fillNoiseCovarianceFlat(R);
        assertEquals(3.5, R[0], TOL);
    }

    public void testFillNoiseCovariance2D() {
        MatrixParameter H = identity2();
        MatrixParameter R = makeMatrix("R", new double[][]{{2.0, 0.5}, {0.5, 3.0}});
        MatrixParameter Y = makeMatrix("Y", new double[][]{{0}, {0}});
        LinearGaussianObservationModel model = new LinearGaussianObservationModel("obs", 2, H, R, Y);

        double[] out = new double[4];
        model.fillNoiseCovarianceFlat(out);
        assertEquals(2.0, out[0], TOL);
        assertEquals(0.5, out[1], TOL);
        assertEquals(0.5, out[2], TOL);
        assertEquals(3.0, out[3], TOL);
    }

    // -------------------------------------------------------------------------
    // fillObservationVector
    // -------------------------------------------------------------------------

    public void testFillObservationVector() {
        LinearGaussianObservationModel model = makeScalar1D(1.0, 1.5, 2.5, 3.5);
        double[] obs = new double[1];

        model.fillObservationVector(0, obs);
        assertEquals(1.5, obs[0], TOL);

        model.fillObservationVector(1, obs);
        assertEquals(2.5, obs[0], TOL);

        model.fillObservationVector(2, obs);
        assertEquals(3.5, obs[0], TOL);
    }

    public void testFillObservationVector2D() {
        MatrixParameter H = identity2();
        MatrixParameter R = identity2();
        // t=0: (10, 20), t=1: (30, 40)
        MatrixParameter Y = makeMatrix("Y", new double[][]{
                {10.0, 30.0},
                {20.0, 40.0}
        });
        LinearGaussianObservationModel model = new LinearGaussianObservationModel("obs", 2, H, R, Y);

        double[] obs = new double[2];
        model.fillObservationVector(0, obs);
        assertEquals(10.0, obs[0], TOL);
        assertEquals(20.0, obs[1], TOL);

        model.fillObservationVector(1, obs);
        assertEquals(30.0, obs[0], TOL);
        assertEquals(40.0, obs[1], TOL);
    }

    public void testFillObservationVectorOutOfBoundsThrows() {
        LinearGaussianObservationModel model = makeScalar1D(1.0, 0.0);
        try {
            model.fillObservationVector(1, new double[1]);
            fail("Expected IllegalArgumentException for out-of-bounds time index");
        } catch (IllegalArgumentException e) {
            // expected
        }
    }

    public void testFillFullyMissingColumnCopiesNaN() {
        MatrixParameter Y = makeMatrix("Y", new double[][]{{Double.NaN}, {Double.NaN}});
        LinearGaussianObservationModel model =
                new LinearGaussianObservationModel("obs", 2, identity2(), identity2(), Y);
        double[] obs = new double[]{1.0, 1.0};
        model.fillObservationVector(0, obs);
        assertTrue(Double.isNaN(obs[0]));
        assertTrue(Double.isNaN(obs[1]));
    }

    // -------------------------------------------------------------------------
    // isObservationColumnMissing
    // -------------------------------------------------------------------------

    public void testIsObservationColumnMissingFalse() {
        LinearGaussianObservationModel model = makeScalar1D(1.0, 5.0);
        assertFalse(model.isObservationColumnMissing(0));
    }

    public void testIsObservationColumnMissingTrue() {
        LinearGaussianObservationModel model = makeScalar1D(1.0, Double.NaN);
        assertTrue(model.isObservationColumnMissing(0));
    }

    public void testMixedMissingAndPresent() {
        LinearGaussianObservationModel model = makeScalar1D(1.0, 1.0, Double.NaN, 3.0);
        assertFalse(model.isObservationColumnMissing(0));
        assertTrue(model.isObservationColumnMissing(1));
        assertFalse(model.isObservationColumnMissing(2));
    }

    public void testIsObservationColumnMissing2D() {
        // t=0 full, t=1 all NaN
        MatrixParameter Y = makeMatrix("Y", new double[][]{
                {1.0, Double.NaN},
                {2.0, Double.NaN}
        });
        LinearGaussianObservationModel model =
                new LinearGaussianObservationModel("obs", 2, identity2(), identity2(), Y);
        assertFalse(model.isObservationColumnMissing(0));
        assertTrue(model.isObservationColumnMissing(1));
    }

    public void testIsObservationColumnMissingOutOfBoundsThrows() {
        LinearGaussianObservationModel model = makeScalar1D(1.0, 0.0);
        try {
            model.isObservationColumnMissing(-1);
            fail("Expected IllegalArgumentException for out-of-bounds time index");
        } catch (IllegalArgumentException e) {
            // expected
        }
    }

    // -------------------------------------------------------------------------
    // Partially observed columns: accepted by the model, rejected by the
    // complete-column accessor used by the moment-form debug engines.
    // -------------------------------------------------------------------------

    private static LinearGaussianObservationModel makePartialColumnModel() {
        MatrixParameter Y = makeMatrix("Y", new double[][]{
                {1.0},
                {Double.NaN}  // only the first coordinate observed at t=0
        });
        return new LinearGaussianObservationModel("obs", 2, identity2(), identity2(), Y);
    }

    public void testPartialColumnIsConstructibleAndNotReportedAsMissing() {
        LinearGaussianObservationModel model = makePartialColumnModel();
        assertEquals(1, model.getTimeCount());
        assertFalse(model.isObservationColumnMissing(0));
        assertEquals(1.0, model.getObservations().getParameterValue(0, 0), TOL);
        assertTrue(Double.isNaN(model.getObservations().getParameterValue(1, 0)));
    }

    public void testPartialMissingFillThrowsAndPointsToCanonicalModes() {
        LinearGaussianObservationModel model = makePartialColumnModel();
        try {
            model.fillObservationVector(0, new double[2]);
            fail("Expected IllegalArgumentException for a partially observed column");
        } catch (IllegalArgumentException e) {
            String message = e.getMessage();
            assertNotNull(message);
            assertTrue(message, message.contains("partially observed"));
            assertTrue(message, message.contains("1 of 2"));
            assertTrue(message, message.contains("canonical"));
            assertTrue(message, message.contains("moment-form"));
        }
    }

    // -------------------------------------------------------------------------
    // Shape validation at construction
    // -------------------------------------------------------------------------

    public void testDesignMatrixRowMismatchThrows() {
        MatrixParameter H = makeMatrix("H", new double[][]{{1.0, 0.0}});   // 1 x 2, need 2 rows
        MatrixParameter Y = makeMatrix("Y", new double[][]{{0}, {0}});
        IllegalArgumentException e = expectConstructionFailure(2, H, identity2(), Y);
        assertTrue(e.getMessage(), e.getMessage().startsWith("designMatrix"));
        assertTrue(e.getMessage(), e.getMessage().contains("1 x 2"));
    }

    public void testNoiseCovarianceNotSquareThrows() {
        MatrixParameter R = makeMatrix("R", new double[][]{{1.0}, {1.0}});  // 2 x 1
        MatrixParameter Y = makeMatrix("Y", new double[][]{{0}, {0}});
        IllegalArgumentException e = expectConstructionFailure(2, identity2(), R, Y);
        assertTrue(e.getMessage(), e.getMessage().startsWith("noiseCovariance"));
        assertTrue(e.getMessage(), e.getMessage().contains("2 x 1"));
    }

    public void testNoiseCovarianceWrongDimensionThrows() {
        MatrixParameter R = makeMatrix("R", new double[][]{{1.0}});  // 1 x 1 for a 2-dimensional observation
        MatrixParameter Y = makeMatrix("Y", new double[][]{{0}, {0}});
        IllegalArgumentException e = expectConstructionFailure(2, identity2(), R, Y);
        assertTrue(e.getMessage(), e.getMessage().startsWith("noiseCovariance"));
    }

    public void testTransposedObservationsThrowWithHint() {
        // Three time points of a 2-dimensional observation written as 3 rows x 2 columns.
        MatrixParameter transposedY = makeMatrix("Y", new double[][]{
                {1.0, 2.0},
                {3.0, 4.0},
                {5.0, 6.0}
        });
        IllegalArgumentException e = expectConstructionFailure(2, identity2(), identity2(), transposedY);
        assertTrue(e.getMessage(), e.getMessage().startsWith("observations"));
        assertTrue(e.getMessage(), e.getMessage().contains("3 x 2"));
        assertTrue(e.getMessage(), e.getMessage().contains("transposed"));
    }

    public void testObservationsRowMismatchThrowsWithoutTransposeHint() {
        MatrixParameter Y = makeMatrix("Y", new double[][]{
                {1.0, 2.0, 3.0, 4.0},
                {1.0, 2.0, 3.0, 4.0},
                {1.0, 2.0, 3.0, 4.0}
        });   // 3 x 4 for a 2-dimensional observation: wrong, but not a transpose
        IllegalArgumentException e = expectConstructionFailure(2, identity2(), identity2(), Y);
        assertTrue(e.getMessage(), e.getMessage().startsWith("observations"));
        assertTrue(e.getMessage(), e.getMessage().contains("3 x 4"));
        assertFalse(e.getMessage(), e.getMessage().contains("transposed"));
    }

    public void testRaggedObservationsThrow() {
        // First column has the right length, second does not: total entries 5 != 2 x 2.
        MatrixParameter ragged = new MatrixParameter("Y", new Parameter[]{
                new Parameter.Default(2),
                new Parameter.Default(3)
        });
        IllegalArgumentException e = expectConstructionFailure(2, identity2(), identity2(), ragged);
        assertTrue(e.getMessage(), e.getMessage().contains("rectangular"));
    }

    public void testObservationsWithoutColumnsThrow() {
        MatrixParameter empty = new MatrixParameter("Y");
        IllegalArgumentException e = expectConstructionFailure(2, identity2(), identity2(), empty);
        assertTrue(e.getMessage(), e.getMessage().contains("at least one time point"));
    }

    // -------------------------------------------------------------------------
    // validateStateDimension
    // -------------------------------------------------------------------------

    public void testValidateStateDimensionAcceptsMatchingColumns() {
        MatrixParameter H = makeMatrix("H", new double[][]{{1.0, 0.0, 0.5}, {0.0, 1.0, -0.5}});
        MatrixParameter Y = makeMatrix("Y", new double[][]{{0}, {0}});
        LinearGaussianObservationModel model = new LinearGaussianObservationModel("obs", 2, H, identity2(), Y);
        model.validateStateDimension(3);
    }

    public void testValidateStateDimensionRejectsMismatch() {
        MatrixParameter H = makeMatrix("H", new double[][]{{1.0, 0.0, 0.5}, {0.0, 1.0, -0.5}});
        MatrixParameter Y = makeMatrix("Y", new double[][]{{0}, {0}});
        LinearGaussianObservationModel model = new LinearGaussianObservationModel("obs", 2, H, identity2(), Y);
        try {
            model.validateStateDimension(2);
            fail("Expected IllegalArgumentException for a 2 x 3 design matrix with a 2-dimensional state");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("2 x 2"));
            assertTrue(e.getMessage(), e.getMessage().contains("2 x 3"));
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
            new LinearGaussianObservationModel("bad", 0, H, R, Y);
            fail("Expected IllegalArgumentException for observationDimension = 0");
        } catch (IllegalArgumentException e) {
            // expected
        }
    }

    public void testNullDesignMatrixThrows() {
        try {
            MatrixParameter R = makeMatrix("R", new double[][]{{1}});
            MatrixParameter Y = makeMatrix("Y", new double[][]{{0}});
            new LinearGaussianObservationModel("bad", 1, null, R, Y);
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
