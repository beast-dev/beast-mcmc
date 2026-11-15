package test.dr.inference.hmc;

import dr.inference.hmc.DiagonalPriorPreconditioner;
import dr.inference.model.PriorPreconditioningProvider;
import dr.xml.Reportable;
import test.dr.math.MathTestCase;

public class DiagonalPriorPreconditionerTest extends MathTestCase {

    public void testReportsFixedStandardDeviations() {
        final DiagonalPriorPreconditioner preconditioner =
                new DiagonalPriorPreconditioner("fixed", new double[]{2.0, 0.5, 1.25});

        assertTrue(preconditioner instanceof PriorPreconditioningProvider);
        assertTrue(preconditioner instanceof Reportable);
        assertEquals(3, preconditioner.getDimension());
        assertEquals(2.0, preconditioner.getStandardDeviation(0), 0.0);
        assertEquals(0.5, preconditioner.getStandardDeviation(1), 0.0);
        assertEquals(1.25, preconditioner.getStandardDeviation(2), 0.0);
        assertTrue(preconditioner.getReport().contains("dimension: 3"));
    }

    public void testCopiesInputArray() {
        final double[] sd = new double[]{1.0, 2.0};
        final DiagonalPriorPreconditioner preconditioner =
                new DiagonalPriorPreconditioner("fixed", sd);
        sd[0] = 10.0;

        assertEquals(1.0, preconditioner.getStandardDeviation(0), 0.0);
    }

    public void testRejectsInvalidStandardDeviations() {
        assertRejects(new double[]{});
        assertRejects(new double[]{1.0, 0.0});
        assertRejects(new double[]{1.0, -1.0});
        assertRejects(new double[]{1.0, Double.NaN});
        assertRejects(new double[]{1.0, Double.POSITIVE_INFINITY});
    }

    private static void assertRejects(final double[] sd) {
        try {
            new DiagonalPriorPreconditioner("bad", sd);
            fail("Expected IllegalArgumentException");
        } catch (IllegalArgumentException expected) {
            // expected
        }
    }
}
