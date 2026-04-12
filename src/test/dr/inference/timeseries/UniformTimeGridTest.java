package test.dr.inference.timeseries;

import dr.inference.timeseries.core.UniformTimeGrid;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Unit tests for {@link UniformTimeGrid}.
 */
public class UniformTimeGridTest extends TestCase {

    public UniformTimeGridTest(String name) {
        super(name);
    }

    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------

    public void testConstructionBasic() {
        UniformTimeGrid grid = new UniformTimeGrid(5, 0.0, 1.0);
        assertEquals(5, grid.getTimeCount());
        assertEquals(0.0, grid.getStartTime(), 0.0);
        assertEquals(1.0, grid.getTimeStep(), 0.0);
    }

    public void testConstructionNonZeroStart() {
        UniformTimeGrid grid = new UniformTimeGrid(3, 2.5, 0.5);
        assertEquals(3, grid.getTimeCount());
        assertEquals(2.5, grid.getStartTime(), 0.0);
        assertEquals(0.5, grid.getTimeStep(), 0.0);
    }

    public void testSingleTimePoint() {
        UniformTimeGrid grid = new UniformTimeGrid(1, 0.0, 1.0);
        assertEquals(1, grid.getTimeCount());
        assertEquals(0.0, grid.getTime(0), 0.0);
    }

    // -------------------------------------------------------------------------
    // getTime(index)
    // -------------------------------------------------------------------------

    public void testGetTimeFirstIndex() {
        UniformTimeGrid grid = new UniformTimeGrid(5, 1.0, 0.5);
        assertEquals(1.0, grid.getTime(0), 1e-15);
    }

    public void testGetTimeLastIndex() {
        UniformTimeGrid grid = new UniformTimeGrid(5, 1.0, 0.5);
        // t[4] = 1.0 + 4 * 0.5 = 3.0
        assertEquals(3.0, grid.getTime(4), 1e-14);
    }

    public void testGetTimeAllIndices() {
        double start = 0.0;
        double step = 0.25;
        int count = 8;
        UniformTimeGrid grid = new UniformTimeGrid(count, start, step);
        for (int i = 0; i < count; i++) {
            double expected = start + i * step;
            assertEquals("time at index " + i, expected, grid.getTime(i), 1e-14);
        }
    }

    // -------------------------------------------------------------------------
    // getDelta(fromIndex, toIndex)
    // -------------------------------------------------------------------------

    public void testGetDeltaAdjacentSteps() {
        UniformTimeGrid grid = new UniformTimeGrid(5, 0.0, 1.0);
        assertEquals(1.0, grid.getDelta(0, 1), 1e-15);
        assertEquals(1.0, grid.getDelta(1, 2), 1e-15);
        assertEquals(1.0, grid.getDelta(3, 4), 1e-15);
    }

    public void testGetDeltaMultiStep() {
        UniformTimeGrid grid = new UniformTimeGrid(6, 0.0, 0.5);
        // delta(0,3) = t[3] - t[0] = 1.5 - 0.0 = 1.5
        assertEquals(1.5, grid.getDelta(0, 3), 1e-14);
    }

    public void testGetDeltaFullRange() {
        UniformTimeGrid grid = new UniformTimeGrid(5, 0.0, 2.0);
        // delta(0,4) = t[4] - t[0] = 8.0
        assertEquals(8.0, grid.getDelta(0, 4), 1e-14);
    }

    // -------------------------------------------------------------------------
    // isRegular()
    // -------------------------------------------------------------------------

    public void testIsRegular() {
        UniformTimeGrid grid = new UniformTimeGrid(10, 0.0, 0.1);
        assertTrue(grid.isRegular());
    }

    // -------------------------------------------------------------------------
    // Invalid construction
    // -------------------------------------------------------------------------

    public void testZeroTimeCountThrows() {
        try {
            new UniformTimeGrid(0, 0.0, 1.0);
            fail("Expected IllegalArgumentException for timeCount = 0");
        } catch (IllegalArgumentException e) {
            // expected
        }
    }

    public void testNegativeTimeCountThrows() {
        try {
            new UniformTimeGrid(-1, 0.0, 1.0);
            fail("Expected IllegalArgumentException for negative timeCount");
        } catch (IllegalArgumentException e) {
            // expected
        }
    }

    public void testZeroTimeStepThrows() {
        try {
            new UniformTimeGrid(5, 0.0, 0.0);
            fail("Expected IllegalArgumentException for timeStep = 0");
        } catch (IllegalArgumentException e) {
            // expected
        }
    }

    public void testNegativeTimeStepThrows() {
        try {
            new UniformTimeGrid(5, 0.0, -0.5);
            fail("Expected IllegalArgumentException for negative timeStep");
        } catch (IllegalArgumentException e) {
            // expected
        }
    }

    // -------------------------------------------------------------------------
    // JUnit 3 boilerplate
    // -------------------------------------------------------------------------

    public static Test suite() {
        return new TestSuite(UniformTimeGridTest.class);
    }
}
