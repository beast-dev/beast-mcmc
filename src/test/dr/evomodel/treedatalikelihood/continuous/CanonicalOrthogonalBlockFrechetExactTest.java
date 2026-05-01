package test.dr.evomodel.treedatalikelihood.continuous;

import dr.inference.model.GivensRotationMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.OrthogonalBlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.Parameter;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockDiagonalSelectionMatrixParameterization;

public class CanonicalOrthogonalBlockFrechetExactTest extends ContinuousTraitTest {

    private static final String FORCE_DENSE_PROPERTY = "beast.experimental.nativeForceDenseAdjointExp";
    private static final double TOL = 1.0e-8;

    public CanonicalOrthogonalBlockFrechetExactTest(final String name) {
        super(name);
    }

    public void testExactFrechetTransitionGradientMatchesDenseFallback() {
        final OrthogonalBlockDiagonalSelectionMatrixParameterization parameterization = buildParameterization();

        final double[] dLogL_dF = new double[]{
                 0.10, -0.25,  0.08,  0.03,
                -0.07,  0.18, -0.11,  0.04,
                 0.06, -0.02,  0.21, -0.09,
                -0.04,  0.05, -0.13,  0.17
        };
        final double[] dLogL_df = new double[]{0.14, -0.09, 0.07, -0.03};
        final double[] stationaryMean = new double[]{0.3, -0.4, 0.8, -0.2};
        final double dt = 0.37;

        final Result exact = computeTransitionGradient(parameterization, stationaryMean, dLogL_dF, dLogL_df, dt, false);
        final Result dense = computeTransitionGradient(parameterization, stationaryMean, dLogL_dF, dLogL_df, dt, true);

        assertVectorEquals("compressed transition gradient", dense.compressed, exact.compressed, TOL);
        assertVectorEquals("rotation transition gradient", dense.rotation, exact.rotation, TOL);
    }

    public void testExactFrechetCovarianceGradientMatchesDenseFallback() {
        final OrthogonalBlockDiagonalSelectionMatrixParameterization parameterization = buildParameterization();
        final MatrixParameter diffusionMatrix = buildDiffusionMatrix();
        final double[] dLogL_dV = new double[]{
                0.40, -0.06,  0.03,  0.01,
               -0.06,  0.28, -0.04,  0.02,
                0.03, -0.04,  0.35, -0.07,
                0.01,  0.02, -0.07,  0.31
        };
        final double dt = 0.29;

        final Result exact = computeCovarianceGradient(parameterization, diffusionMatrix, dLogL_dV, dt, false);
        final Result dense = computeCovarianceGradient(parameterization, diffusionMatrix, dLogL_dV, dt, true);

        assertVectorEquals("compressed covariance gradient", dense.compressed, exact.compressed, TOL);
        assertVectorEquals("rotation covariance gradient", dense.rotation, exact.rotation, TOL);
    }

    public void testExactFrechetNearNilpotentTransitionGradientMatchesDenseFallback() {
        final double theta1 = 0.60;
        final double theta2 = -0.45;
        final double rho1 = 1.0;
        final double rho2 = 0.9;
        final double small = 1.0e-8;

        final OrthogonalBlockDiagonalSelectionMatrixParameterization parameterization = buildParameterization(
                new double[]{rho1, rho2},
                new double[]{theta1, theta2},
                new double[]{
                        rho1 * Math.sin(theta1) - small,
                        -rho2 * Math.sin(theta2) + small
                });

        final double[] dLogL_dF = new double[]{
                 0.12, -0.09,  0.04,  0.01,
                -0.05,  0.17, -0.08,  0.03,
                 0.03, -0.07,  0.19, -0.06,
                -0.02,  0.04, -0.05,  0.14
        };
        final double[] dLogL_df = new double[]{0.11, -0.07, 0.05, -0.02};
        final double[] stationaryMean = new double[]{0.2, -0.3, 0.6, -0.1};
        final double dt = 0.41;

        final Result exact = computeTransitionGradient(parameterization, stationaryMean, dLogL_dF, dLogL_df, dt, false);
        final Result dense = computeTransitionGradient(parameterization, stationaryMean, dLogL_dF, dLogL_df, dt, true);

        assertVectorEquals("near-nilpotent compressed transition gradient", dense.compressed, exact.compressed, TOL);
        assertVectorEquals("near-nilpotent rotation transition gradient", dense.rotation, exact.rotation, TOL);
    }

    private Result computeTransitionGradient(final OrthogonalBlockDiagonalSelectionMatrixParameterization parameterization,
                                             final double[] stationaryMean,
                                             final double[] dLogL_dF,
                                             final double[] dLogL_df,
                                             final double dt,
                                             final boolean forceDense) {
        final String previous = System.getProperty(FORCE_DENSE_PROPERTY);
        try {
            if (forceDense) {
                System.setProperty(FORCE_DENSE_PROPERTY, "true");
            } else {
                System.clearProperty(FORCE_DENSE_PROPERTY);
            }

            final Result result = new Result(10, 4);
            parameterization.accumulateNativeGradientFromTransitionFlat(
                    dt, stationaryMean, dLogL_dF, dLogL_df, result.compressed, result.rotation);
            return result;
        } finally {
            restoreProperty(previous);
        }
    }

    private Result computeCovarianceGradient(final OrthogonalBlockDiagonalSelectionMatrixParameterization parameterization,
                                             final MatrixParameter diffusionMatrix,
                                             final double[] dLogL_dV,
                                             final double dt,
                                             final boolean forceDense) {
        final String previous = System.getProperty(FORCE_DENSE_PROPERTY);
        try {
            if (forceDense) {
                System.setProperty(FORCE_DENSE_PROPERTY, "true");
            } else {
                System.clearProperty(FORCE_DENSE_PROPERTY);
            }

            final Result result = new Result(10, 4);
            parameterization.accumulateNativeGradientFromCovarianceStationaryFlat(
                    diffusionMatrix, dt, dLogL_dV, result.compressed, result.rotation);
            return result;
        } finally {
            restoreProperty(previous);
        }
    }

    private void restoreProperty(final String previous) {
        if (previous == null) {
            System.clearProperty(FORCE_DENSE_PROPERTY);
        } else {
            System.setProperty(FORCE_DENSE_PROPERTY, previous);
        }
    }

    private OrthogonalBlockDiagonalSelectionMatrixParameterization buildParameterization() {
        return buildParameterization(
                new double[]{0.85, 1.10},
                new double[]{0.55, -0.35},
                new double[]{0.27, -0.19});
    }

    private OrthogonalBlockDiagonalSelectionMatrixParameterization buildParameterization(final double[] rhoValues,
                                                                                         final double[] thetaValues,
                                                                                         final double[] tValues) {
        final Parameter angles = new Parameter.Default("angles.frechetExact",
                new double[]{0.10, -0.12, 0.08, 0.04, -0.06, 0.11});
        final GivensRotationMatrixParameter rotation =
                new GivensRotationMatrixParameter("rotation.frechetExact", 4, angles);
        final Parameter scalar = new Parameter.Default("scalar.frechetExact", new double[0]);
        final Parameter rho = new Parameter.Default("rho.frechetExact", rhoValues);
        final Parameter theta = new Parameter.Default("theta.frechetExact", thetaValues);
        final Parameter t = new Parameter.Default("t.frechetExact", tValues);

        final OrthogonalBlockDiagonalPolarStableMatrixParameter blockSelection =
                new OrthogonalBlockDiagonalPolarStableMatrixParameter(
                        "selection.frechetExact",
                        rotation,
                        scalar,
                        rho,
                        theta,
                        t);

        return new OrthogonalBlockDiagonalSelectionMatrixParameterization(blockSelection, rotation);
    }

    private MatrixParameter buildDiffusionMatrix() {
        final MatrixParameter q = new MatrixParameter("q.frechetExact", 4, 4);
        final double[][] values = new double[][]{
                {2.0, 0.20, 0.05, 0.01},
                {0.20, 1.7, 0.12, -0.03},
                {0.05, 0.12, 1.9, 0.18},
                {0.01, -0.03, 0.18, 1.6}
        };
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                q.setParameterValueQuietly(i, j, values[i][j]);
            }
        }
        q.fireParameterChangedEvent();
        return q;
    }

    private void assertVectorEquals(final String label,
                                    final double[] expected,
                                    final double[] actual,
                                    final double tolerance) {
        assertEquals(label + " length", expected.length, actual.length);
        for (int i = 0; i < expected.length; ++i) {
            assertEquals(label + "[" + i + "]", expected[i], actual[i], tolerance);
        }
    }

    private void assertMatrixEquals(final String label,
                                    final double[][] expected,
                                    final double[][] actual,
                                    final double tolerance) {
        assertEquals(label + " rows", expected.length, actual.length);
        for (int i = 0; i < expected.length; ++i) {
            assertEquals(label + " cols[" + i + "]", expected[i].length, actual[i].length);
            for (int j = 0; j < expected[i].length; ++j) {
                assertEquals(label + "[" + i + "][" + j + "]", expected[i][j], actual[i][j], tolerance);
            }
        }
    }

    private static final class Result {
        private final double[] compressed;
        private final double[] rotation;

        private Result(final int compressedDimension,
                       final int matrixDimension) {
            this.compressed = new double[compressedDimension];
            this.rotation = new double[matrixDimension * matrixDimension];
        }
    }
}
