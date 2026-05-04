package test.dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalFrechetHelper;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalExpSolver;
import dr.inference.model.GivensRotationMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.OrthogonalBlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.Parameter;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockDiagonalSelectionMatrixParameterization;
import org.ejml.data.DenseMatrix64F;

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

    public void testOrthogonalPolarBlocksUseEqualDiagonalExactPlanKernels() {
        final OrthogonalBlockDiagonalSelectionMatrixParameterization parameterization = buildParameterization();
        final double[] dLogL_dF = new double[]{
                 0.10, -0.25,  0.08,  0.03,
                -0.07,  0.18, -0.11,  0.04,
                 0.06, -0.02,  0.21, -0.09,
                -0.04,  0.05, -0.13,  0.17
        };
        final double[] dLogL_df = new double[]{0.14, -0.09, 0.07, -0.03};
        final double[] stationaryMean = new double[]{0.3, -0.4, 0.8, -0.2};

        BlockDiagonalFrechetHelper.resetExactPlanInstrumentation();
        computeTransitionGradient(parameterization, stationaryMean, dLogL_dF, dLogL_df, 0.37, false);

        assertTrue("exact plan updates parameters",
                BlockDiagonalFrechetHelper.getExactPlanParameterUpdateCount() > 0L);
        assertTrue("exact plan evaluates at branch length",
                BlockDiagonalFrechetHelper.getExactPlanTimeEvaluationCount() > 0L);
        assertTrue("orthogonal polar blocks use equal-diagonal kernels",
                BlockDiagonalFrechetHelper.getExactPlanEqualDiagonalEvaluationCount() > 0L);
        assertTrue("positive-root polar blocks use real distinct coefficient fast path",
                BlockDiagonalFrechetHelper.getExactPlanCoefficientDistinctRealEvaluationCount() > 0L);
        assertEquals("all distinct positive-root coefficients use real fast path",
                BlockDiagonalFrechetHelper.getExactPlanCoefficientDistinctEvaluationCount(),
                BlockDiagonalFrechetHelper.getExactPlanCoefficientDistinctRealEvaluationCount());
        assertEquals("orthogonal polar blocks avoid generic Sylvester solves",
                0L,
                BlockDiagonalFrechetHelper.getExactPlanGenericSolve4x4CallCount());
    }

    public void testExactFrechetPlanReusesSameBranchTime() {
        final BlockDiagonalExpSolver.BlockStructure structure =
                new BlockDiagonalExpSolver.BlockStructure(4, new int[]{0, 2}, new int[]{2, 2});
        final BlockDiagonalFrechetHelper helper = new BlockDiagonalFrechetHelper(structure);
        final double[] blockDParams = new double[]{
                0.93, 0.93, 1.04, 1.04,
                0.18, 0.0, -0.14,
                -0.11, 0.0, 0.09
        };
        final DenseMatrix64F firstInput = new DenseMatrix64F(4, 4);
        final DenseMatrix64F secondInput = new DenseMatrix64F(4, 4);
        final DenseMatrix64F firstOut = new DenseMatrix64F(4, 4);
        final DenseMatrix64F secondOut = new DenseMatrix64F(4, 4);
        for (int i = 0; i < firstInput.data.length; ++i) {
            firstInput.data[i] = 0.01 * (i + 1);
            secondInput.data[i] = -0.02 * (i + 1);
        }

        BlockDiagonalFrechetHelper.resetExactPlanInstrumentation();
        helper.frechetAdjointExpInDBasis(blockDParams, firstInput, 0.37, firstOut);
        helper.frechetAdjointExpInDBasis(blockDParams, secondInput, 0.37, secondOut);

        assertEquals("same branch time evaluates coefficients once",
                1L,
                BlockDiagonalFrechetHelper.getExactPlanTimeEvaluationCount());
        assertEquals("second same branch time call hits cache",
                1L,
                BlockDiagonalFrechetHelper.getExactPlanCacheHitCount());
        assertEquals("both inputs still apply the exact plan",
                2L,
                BlockDiagonalFrechetHelper.getExactPlanApplicationCount());
        assertEquals("same branch time fills each 2x2 block-pair coefficient once",
                4L,
                BlockDiagonalFrechetHelper.getExactPlanKernel2x2EqualDiagonalEvaluationCount());
        assertEquals("same branch time does not use generic 2x2 kernels",
                0L,
                BlockDiagonalFrechetHelper.getExactPlanKernel2x2GenericEvaluationCount());
    }

    public void testExactFrechetPlanRoutesComplexDistinctShapes() {
        final BlockDiagonalExpSolver.BlockStructure structure =
                new BlockDiagonalExpSolver.BlockStructure(4, new int[]{0, 2}, new int[]{2, 2});
        final BlockDiagonalFrechetHelper helper = new BlockDiagonalFrechetHelper(structure);
        final double[] blockDParams = new double[]{
                0.93, 0.93, 1.04, 1.04,
                0.55, 0.0, 0.48,
                -0.41, 0.0, 0.36
        };
        final DenseMatrix64F input = new DenseMatrix64F(4, 4);
        final DenseMatrix64F out = new DenseMatrix64F(4, 4);
        for (int i = 0; i < input.data.length; ++i) {
            input.data[i] = 0.03 * (i + 1);
        }

        BlockDiagonalFrechetHelper.resetExactPlanInstrumentation();
        helper.computeForwardFrechetInDBasis(blockDParams, input, 0.37, out);

        assertEquals("one real/real distinct 2x2 block pair",
                1L,
                BlockDiagonalFrechetHelper.getExactPlanCoefficientDistinctRealEvaluationCount());
        assertEquals("one imaginary/real distinct 2x2 block pair",
                1L,
                BlockDiagonalFrechetHelper.getExactPlanCoefficientDistinctLeftImagRightRealEvaluationCount());
        assertEquals("one real/imaginary distinct 2x2 block pair",
                1L,
                BlockDiagonalFrechetHelper.getExactPlanCoefficientDistinctLeftRealRightImagEvaluationCount());
        assertEquals("one imaginary/imaginary distinct 2x2 block pair",
                1L,
                BlockDiagonalFrechetHelper.getExactPlanCoefficientDistinctBothImagEvaluationCount());
        assertEquals("all block pairs used equal-diagonal distinct coefficients",
                4L,
                BlockDiagonalFrechetHelper.getExactPlanCoefficientDistinctEvaluationCount());
        for (int i = 0; i < out.data.length; ++i) {
            assertTrue("finite exact Frechet output[" + i + "]", Double.isFinite(out.data[i]));
        }
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
