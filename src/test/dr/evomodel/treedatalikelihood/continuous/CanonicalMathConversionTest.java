package test.dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.continuous.canonical.math.CanonicalInversionResult;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.GaussianFormConverter;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.inference.model.GivensRotationMatrixParameter;
import dr.inference.model.Parameter;
import junit.framework.TestCase;
import org.ejml.data.DenseMatrix64F;

public final class CanonicalMathConversionTest extends TestCase {

    private static final double TOL = 1.0e-11;

    public void testMomentCanonicalStateRoundTrip() {
        final int dim = 3;
        final double[] mean = {1.25, -0.75, 2.5};
        final double[] covariance = {
                3.0, 0.4, 0.2,
                0.4, 2.5, 0.3,
                0.2, 0.3, 1.75
        };
        final CanonicalGaussianState state = new CanonicalGaussianState(dim);
        final GaussianFormConverter.Workspace workspace = new GaussianFormConverter.Workspace();
        workspace.ensureDim(dim);

        GaussianFormConverter.fillStateFromMoments(mean, covariance, dim, workspace, state);

        final double[] recoveredMean = new double[dim];
        final double[] recoveredCovariance = new double[dim * dim];
        GaussianFormConverter.fillMomentsFromState(
                state, recoveredMean, recoveredCovariance, dim, workspace);

        assertArrayClose("mean round trip", mean, recoveredMean, TOL);
        assertArrayClose("covariance round trip", covariance, recoveredCovariance, TOL);
    }

    public void testFlatAndMatrixTransitionMomentConversionMatch() {
        final int dim = 2;
        final double[] flatTransition = {
                0.9, 0.15,
                -0.05, 0.8
        };
        final double[][] matrixTransition = {
                {0.9, 0.15},
                {-0.05, 0.8}
        };
        final double[] offset = {0.2, -0.35};
        final double[] flatOmega = {
                1.4, 0.25,
                0.25, 1.1
        };
        final double[][] matrixOmega = {
                {1.4, 0.25},
                {0.25, 1.1}
        };
        final CanonicalGaussianTransition flat = new CanonicalGaussianTransition(dim);
        final CanonicalGaussianTransition matrix = new CanonicalGaussianTransition(dim);
        final GaussianFormConverter.Workspace workspace = new GaussianFormConverter.Workspace();
        workspace.ensureDim(dim);

        GaussianFormConverter.fillTransitionFromMoments(
                flatTransition, offset, flatOmega, dim, workspace, flat);
        GaussianFormConverter.fillTransitionFromMoments(
                matrixTransition, offset, matrixOmega, matrix);

        assertArrayClose("precisionXX", flat.precisionXX, matrix.precisionXX, TOL);
        assertArrayClose("precisionXY", flat.precisionXY, matrix.precisionXY, TOL);
        assertArrayClose("precisionYX", flat.precisionYX, matrix.precisionYX, TOL);
        assertArrayClose("precisionYY", flat.precisionYY, matrix.precisionYY, TOL);
        assertArrayClose("informationX", flat.informationX, matrix.informationX, TOL);
        assertArrayClose("informationY", flat.informationY, matrix.informationY, TOL);
        assertEquals("log normalizer", flat.logNormalizer, matrix.logNormalizer, TOL);
    }

    public void testCompactAndGeneralSpdInversionMatch() {
        final int dim = 3;
        final double[] source = {
                4.0, 0.4, 0.2,
                0.4, 3.0, 0.5,
                0.2, 0.5, 2.0
        };
        final double[] compactSource = source.clone();
        final double[] compactInverse = new double[dim * dim];
        final double[] generalInverse = new double[dim * dim];

        final double compactLogDet = MatrixOps.invertSPDCompact(
                compactSource,
                compactInverse,
                dim,
                new double[dim],
                new double[dim]);
        MatrixOps.invertSPD(source.clone(), generalInverse, dim);
        final double generalLogDet = MatrixOps.logDeterminant(source.clone(), dim);

        assertArrayClose("compact vs general inverse", generalInverse, compactInverse, TOL);
        assertEquals("compact vs general log determinant", generalLogDet, compactLogDet, TOL);
    }

    public void testEJMLBoundaryConversionUsesFlatRowMajorStorage() {
        final int dim = 2;
        final DenseMatrix64F matrix = new DenseMatrix64F(dim, dim);
        matrix.set(0, 0, 1.0);
        matrix.set(0, 1, 2.0);
        matrix.set(1, 0, 3.0);
        matrix.set(1, 1, 4.0);
        final double[] flat = new double[dim * dim];

        MatrixOps.toFlat(matrix, flat, dim);
        assertArrayClose("ejml to flat", new double[] {1.0, 2.0, 3.0, 4.0}, flat, TOL);

        final DenseMatrix64F roundTrip = new DenseMatrix64F(dim, dim);
        MatrixOps.fromFlat(flat, roundTrip, dim);
        assertEquals(1.0, roundTrip.get(0, 0), TOL);
        assertEquals(3.0, roundTrip.get(1, 0), TOL);
        assertEquals(4.0, roundTrip.get(1, 1), TOL);
    }

    public void testGivensFlatPullbackCanFillExistingBuffer() {
        final int dim = 3;
        final Parameter angles = new Parameter.Default(new double[] {0.15, -0.30, 0.45});
        final GivensRotationMatrixParameter rotation =
                new GivensRotationMatrixParameter("rotation.fillPullback", dim, angles);
        final double[] gradient = {
                0.4, -0.2, 0.7,
                0.1, 0.5, -0.3,
                -0.6, 0.8, 0.2
        };
        final double[] expected = new double[angles.getDimension()];
        final double[] actual = {-99.0, -99.0, -99.0, -99.0, -99.0};

        rotation.fillPullBackGradientFlat(gradient, dim, expected);
        rotation.fillPullBackGradientFlat(gradient, dim, actual, 1);

        assertEquals(-99.0, actual[0], TOL);
        assertEquals(-99.0, actual[4], TOL);
        for (int i = 0; i < expected.length; ++i) {
            assertEquals("angle gradient " + i, expected[i], actual[i + 1], TOL);
        }
    }

    public void testSafeInvertPrecisionHandlesMissingDiagonalLocally() {
        final int dim = 3;
        final double[] precision = {
                2.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                0.0, 0.0, 4.0
        };
        final double[] variance = new double[dim * dim];

        final CanonicalInversionResult result = MatrixOps.safeInvertPrecision(precision, variance, dim);

        assertEquals(CanonicalInversionResult.Code.PARTIALLY_OBSERVED, result.getReturnCode());
        assertEquals(2, result.getEffectiveDimension());
        assertEquals(Math.log(8.0), result.getLogDeterminant(), TOL);
        assertEquals(0.5, variance[0], TOL);
        assertEquals(Double.POSITIVE_INFINITY, variance[4]);
        assertEquals(0.25, variance[8], TOL);
        assertEquals(0.0, variance[1], TOL);
        assertEquals(0.0, variance[5], TOL);
    }

    public void testSafeInvertVarianceHandlesMissingDiagonalLocally() {
        final int dim = 3;
        final double[] variance = {
                2.0, 0.0, 0.0,
                0.0, Double.POSITIVE_INFINITY, 0.0,
                0.0, 0.0, 4.0
        };
        final double[] precision = new double[dim * dim];

        final CanonicalInversionResult result = MatrixOps.safeInvertVariance(variance, precision, dim);

        assertEquals(CanonicalInversionResult.Code.PARTIALLY_OBSERVED, result.getReturnCode());
        assertEquals(2, result.getEffectiveDimension());
        assertEquals(Math.log(8.0), result.getLogDeterminant(), TOL);
        assertEquals(0.5, precision[0], TOL);
        assertEquals(0.0, precision[4], TOL);
        assertEquals(0.25, precision[8], TOL);
        assertEquals(0.0, precision[1], TOL);
        assertEquals(0.0, precision[5], TOL);
    }

    public void testSafeInvertVariancePreservesCollapsedDiagonal() {
        final int dim = 2;
        final double[] variance = {
                0.0, 0.0,
                0.0, Double.POSITIVE_INFINITY
        };
        final double[] precision = new double[dim * dim];

        final CanonicalInversionResult result = MatrixOps.safeInvertVariance(variance, precision, dim);

        assertEquals(CanonicalInversionResult.Code.PARTIALLY_OBSERVED, result.getReturnCode());
        assertEquals(1, result.getEffectiveDimension());
        assertEquals(Double.POSITIVE_INFINITY, result.getLogDeterminant());
        assertEquals(Double.POSITIVE_INFINITY, precision[0]);
        assertEquals(0.0, precision[3], TOL);
    }

    public void testMatrixOpsMatMulMatchesReferenceAcrossDimensions() {
        assertMatMulMatchesReference(3);
        assertMatMulMatchesReference(7);
        assertMatMulMatchesReference(8);
    }

    public void testMatrixOpsMatMulTransposedRightMatchesReferenceAcrossDimensions() {
        assertMatMulTransposedRightMatchesReference(3);
        assertMatMulTransposedRightMatchesReference(7);
        assertMatMulTransposedRightMatchesReference(8);
    }

    private void assertMatMulMatchesReference(final int dim) {
        final double[] left = fillPattern(dim, 0.17, -0.03);
        final double[] right = fillPattern(dim, -0.11, 0.05);
        final double[] expected = new double[dim * dim];
        final double[] actual = new double[dim * dim];

        referenceMatMul(left, right, expected, dim);
        MatrixOps.matMul(left, right, actual, dim);

        assertArrayClose("matMul dim=" + dim, expected, actual, TOL);
    }

    private void assertMatMulTransposedRightMatchesReference(final int dim) {
        final double[] left = fillPattern(dim, 0.13, 0.07);
        final double[] right = fillPattern(dim, -0.19, 0.02);
        final double[] expected = new double[dim * dim];
        final double[] actual = new double[dim * dim];

        referenceMatMulTransposedRight(left, right, expected, dim);
        MatrixOps.matMulTransposedRight(left, right, actual, dim);

        assertArrayClose("matMulTransposedRight dim=" + dim, expected, actual, TOL);
    }

    public void testSymmetricSandwichTransposeLeftMatchesReference() {
        assertSymmetricSandwichTransposeLeftMatchesReference(4);
        assertSymmetricSandwichTransposeLeftMatchesReference(8);
    }

    public void testSymmetricSandwichTransposeRightMatchesReference() {
        assertSymmetricSandwichTransposeRightMatchesReference(4);
        assertSymmetricSandwichTransposeRightMatchesReference(8);
    }

    private void assertSymmetricSandwichTransposeLeftMatchesReference(final int dim) {
        final double[] left = fillPattern(dim, 0.09, -0.04);
        final double[] middle = fillSymmetricPattern(dim);
        final double[] middleTimesLeft = new double[dim * dim];
        final double[] expectedScratch = new double[dim * dim];
        final double[] expected = new double[dim * dim];
        final double[] actual = new double[dim * dim];

        referenceMatMul(middle, left, expectedScratch, dim);
        referenceMatMulTransposedLeft(left, expectedScratch, expected, dim);
        MatrixOps.symmetricSandwichTransposeLeft(left, middle, actual, middleTimesLeft, dim);

        assertArrayClose("symmetric transpose-left sandwich scratch dim=" + dim,
                expectedScratch, middleTimesLeft, TOL);
        assertArrayClose("symmetric transpose-left sandwich dim=" + dim, expected, actual, TOL);
    }

    private void assertSymmetricSandwichTransposeRightMatchesReference(final int dim) {
        final double[] left = fillPattern(dim, -0.07, 0.06);
        final double[] middle = fillSymmetricPattern(dim);
        final double[] leftTimesMiddle = new double[dim * dim];
        final double[] expectedScratch = new double[dim * dim];
        final double[] expected = new double[dim * dim];
        final double[] actual = new double[dim * dim];

        referenceMatMul(left, middle, expectedScratch, dim);
        referenceMatMulTransposedRight(expectedScratch, left, expected, dim);
        MatrixOps.symmetricSandwichTransposeRight(left, middle, actual, leftTimesMiddle, dim);

        assertArrayClose("symmetric transpose-right sandwich scratch dim=" + dim,
                expectedScratch, leftTimesMiddle, TOL);
        assertArrayClose("symmetric transpose-right sandwich dim=" + dim, expected, actual, TOL);
    }

    private static void assertArrayClose(final String label,
                                         final double[] expected,
                                         final double[] actual,
                                         final double tolerance) {
        assertEquals(label + " length", expected.length, actual.length);
        for (int i = 0; i < expected.length; ++i) {
            assertEquals(label + "[" + i + "]", expected[i], actual[i], tolerance);
        }
    }

    private static double[] fillPattern(final int dim,
                                        final double rowScale,
                                        final double columnScale) {
        final double[] values = new double[dim * dim];
        for (int i = 0; i < dim; ++i) {
            final int rowOffset = i * dim;
            for (int j = 0; j < dim; ++j) {
                values[rowOffset + j] =
                        rowScale * (i + 1) + columnScale * (j + 1) + 0.01 * ((i + j) % 3);
            }
        }
        return values;
    }

    private static double[] fillSymmetricPattern(final int dim) {
        final double[] values = new double[dim * dim];
        for (int i = 0; i < dim; ++i) {
            for (int j = i; j < dim; ++j) {
                final double value = (i == j ? 2.0 : 0.0) + 0.03 * (i + 1) - 0.02 * (j + 1);
                values[i * dim + j] = value;
                values[j * dim + i] = value;
            }
        }
        return values;
    }

    private static void referenceMatMul(final double[] left,
                                        final double[] right,
                                        final double[] out,
                                        final int dim) {
        for (int i = 0; i < dim; ++i) {
            final int rowOffset = i * dim;
            for (int j = 0; j < dim; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dim; ++k) {
                    sum += left[rowOffset + k] * right[k * dim + j];
                }
                out[rowOffset + j] = sum;
            }
        }
    }

    private static void referenceMatMulTransposedRight(final double[] left,
                                                       final double[] right,
                                                       final double[] out,
                                                       final int dim) {
        for (int i = 0; i < dim; ++i) {
            final int rowOffset = i * dim;
            for (int j = 0; j < dim; ++j) {
                double sum = 0.0;
                final int rightOffset = j * dim;
                for (int k = 0; k < dim; ++k) {
                    sum += left[rowOffset + k] * right[rightOffset + k];
                }
                out[rowOffset + j] = sum;
            }
        }
    }

    private static void referenceMatMulTransposedLeft(final double[] left,
                                                      final double[] right,
                                                      final double[] out,
                                                      final int dim) {
        for (int i = 0; i < dim; ++i) {
            final int rowOffset = i * dim;
            for (int j = 0; j < dim; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dim; ++k) {
                    sum += left[k * dim + i] * right[k * dim + j];
                }
                out[rowOffset + j] = sum;
            }
        }
    }
}
