package test.dr.inference.model;

import dr.inference.model.AbstractComputedCompoundMatrix;
import dr.inference.model.Parameter;
import test.dr.math.MathTestCase;

import java.util.Collections;

/**
 * Regression for stale computed-matrix caches after reject/restore.
 *
 * The computed matrix cache must be invalidated on restore because the cache
 * buffer is owned by subclasses and cannot be restored generically in the
 * abstract base class.
 */
public class AbstractComputedCompoundMatrixRestoreTest extends MathTestCase {

    private static final class OneByOneComputedMatrix extends AbstractComputedCompoundMatrix {

        private final Parameter source;
        private final double[] cached = new double[1];

        private OneByOneComputedMatrix(final Parameter source) {
            super("oneByOneComputed", 1, Collections.singletonList(source));
            this.source = source;
        }

        @Override
        protected double computeEntry(final int row, final int col) {
            return source.getParameterValue(0);
        }

        @Override
        protected double getCachedValue(final int row, final int col) {
            return cached[0];
        }

        @Override
        protected void updateCache() {
            cached[0] = computeEntry(0, 0);
        }
    }

    public void testRestoreInvalidatesComputedCache() {

        final Parameter.Default source = new Parameter.Default(1.0);
        final OneByOneComputedMatrix computed = new OneByOneComputedMatrix(source);

        assertTrue(computed.isImmutable());

        // Ensure computed cache is initialized and marked known before store.
        assertEquals(1.0, computed.getParameterValue(0, 0), 1e-15);

        computed.storeParameterValues();
        source.storeParameterValues();

        // Proposed state.
        source.setParameterValue(0, 2.0);
        assertEquals(2.0, computed.getParameterValue(0, 0), 1e-15);

        // Reject proposal: restore child first, then computed wrapper.
        source.restoreParameterValues();
        computed.restoreParameterValues();

        // Must reflect restored child state, not stale proposed cached value.
        assertEquals(1.0, computed.getParameterValue(0, 0), 1e-15);
    }
}
