package test.dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalNumerics;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalNumericsOptions;
import junit.framework.TestCase;
import org.ejml.data.DenseMatrix64F;

public class CanonicalNumericsTest extends TestCase {

    public CanonicalNumericsTest(final String name) {
        super(name);
    }

    public void testStrictDenseSpdFallbackUsesJitterPolicy() {
        final DenseMatrix64F source = new DenseMatrix64F(new double[][] {
                {0.0, 0.0},
                {0.0, 2.0}
        });
        final DenseMatrix64F inverse = new DenseMatrix64F(2, 2);

        CanonicalNumerics.safeInvertSymmetricPositiveDefinite(
                source,
                inverse,
                new DenseMatrix64F(2, 2),
                new DenseMatrix64F(2, 2),
                new double[2][2],
                new double[2][2],
                new double[2][2],
                new double[2][2],
                new CanonicalNumericsOptions(1.0e-6, 1.0e-6, 0, 2, true, false, false),
                "strict-fallback-test",
                null);

        assertTrue(CanonicalNumerics.isFinite(inverse));
        assertEquals(inverse.unsafe_get(0, 1), inverse.unsafe_get(1, 0), 0.0);
        assertTrue("jittered zero pivot should become invertible", inverse.unsafe_get(0, 0) > 1.0e5);
        assertEquals(0.5, inverse.unsafe_get(1, 1), 1.0e-6);
    }

    public void testOrthogonalBlockFlatInversionAppliesPivotFloor() {
        final DenseMatrix64F source = new DenseMatrix64F(new double[][] {
                {0.0, 0.0},
                {0.0, 4.0}
        });
        final double[] covariance = new double[4];
        final double[] inverse = new double[4];

        final double logDet = CanonicalNumerics.copyAndInvertPositiveDefiniteFlat(
                source,
                covariance,
                inverse,
                new double[4],
                new double[4],
                new CanonicalNumericsOptions(1.0e-6, 1.0e-6, 8, 8, false, false, false));

        assertTrue(Double.isFinite(logDet));
        assertTrue(CanonicalNumerics.isFinite(covariance));
        assertTrue(CanonicalNumerics.isFinite(inverse));
        assertEquals(0.0, covariance[0], 0.0);
        assertEquals(4.0, covariance[3], 0.0);
        assertTrue("pivot floor should produce a finite inverse for the clipped pivot", inverse[0] > 1.0e5);
        assertEquals(0.25, inverse[3], 1.0e-9);
    }

    public void testFiniteChecksAndDenseMatrixSummarySeparateNaNAndInfinity() {
        final DenseMatrix64F matrix = new DenseMatrix64F(new double[][] {
                {Double.NaN, Double.POSITIVE_INFINITY},
                {Double.NEGATIVE_INFINITY, 4.0}
        });

        assertFalse(CanonicalNumerics.isFinite(matrix));
        assertTrue(CanonicalNumerics.hasNaN(matrix));
        assertTrue(CanonicalNumerics.hasInfinity(matrix));
        final String summary = CanonicalNumerics.summarizeDenseMatrix(matrix);
        assertTrue(summary.contains("nan=1"));
        assertTrue(summary.contains("posInf=1"));
        assertTrue(summary.contains("negInf=1"));
        assertTrue(summary.contains("minFinite=4.0"));
        assertTrue(summary.contains("maxFinite=4.0"));
    }

    public void testSpdFailureDumpCallbackRunsOnFailedInversion() {
        final DenseMatrix64F source = new DenseMatrix64F(new double[][] {
                {Double.NaN, 0.0},
                {0.0, 1.0}
        });
        final boolean[] callbackInvoked = {false};

        try {
            CanonicalNumerics.safeInvertSymmetricPositiveDefinite(
                    source,
                    new DenseMatrix64F(2, 2),
                    new DenseMatrix64F(2, 2),
                    new DenseMatrix64F(2, 2),
                    new double[2][2],
                    new double[2][2],
                    new double[2][2],
                    new double[2][2],
                    new CanonicalNumericsOptions(1.0e-12, 1.0e-12, 0, 0, false, true, false),
                    "dump-callback-test",
                    (sb, context, originalSource, symmetrizedSource, jitterBase) -> {
                        callbackInvoked[0] = true;
                        sb.append(",\n\"testContext\":\"").append(context).append("\"\n");
                    });
            fail("Expected failed SPD inversion");
        } catch (final IllegalStateException expected) {
            assertTrue("debug dump callback should be invoked before throwing", callbackInvoked[0]);
        }
    }
}
