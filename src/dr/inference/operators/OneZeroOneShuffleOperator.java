package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.math.MathUtils;

/**
 * Moves which TWO labeled indices are the extremes:
 *  - extremeIndex[0] = index of the alpha entry that is forced to be 0
 *  - extremeIndex[1] = index of the alpha entry that is forced to be 1
 *
 * This operator:
 *  - chooses whether to move the 0-label or the 1-label
 *  - chooses a destination index among the NON-extreme indices
 *  - swaps the label (extremeIndex entry) with dst
 *
 */

/*
* @author: Filippo Monti
 */
public final class OneZeroOneShuffleOperator extends SimpleMCMCOperator {

    private final Parameter alphaRates;    // length K
    private final Parameter extremeIndex;  // length 2, integer-like doubles
    private final double tol;

    public OneZeroOneShuffleOperator(final Parameter alphaRates,
                                     final Parameter extremeIndex,
                                     final double weight,
                                     final double tol) {
        if (alphaRates == null) throw new IllegalArgumentException("alphaRates must be non-null");
        if (extremeIndex == null) throw new IllegalArgumentException("extremeIndex must be non-null");
        if (alphaRates.getDimension() < 2) throw new IllegalArgumentException("alphaRates dimension must be >= 2");
        if (extremeIndex.getDimension() != 2) throw new IllegalArgumentException("extremeIndex dimension must be 2");
        if (tol < 0.0) throw new IllegalArgumentException("tol must be >= 0");

        this.alphaRates = alphaRates;
        this.extremeIndex = extremeIndex;
        this.tol = tol;

        setWeight(weight);
    }

    @Override
    public String getOperatorName() {
        return "oneZeroOneShuffleOperator(" + alphaRates.getParameterName() + ")";
    }

    @Override
    public double doOperation() {

        final int K = alphaRates.getDimension();

        final int i0 = toIndexStrict(extremeIndex.getParameterValue(0), "extremeIndex[0]");
        final int i1 = toIndexStrict(extremeIndex.getParameterValue(1), "extremeIndex[1]");

        if (i0 < 0 || i0 >= K) throw new IllegalArgumentException("extremeIndex[0] out of range: " + i0);
        if (i1 < 0 || i1 >= K) throw new IllegalArgumentException("extremeIndex[1] out of range: " + i1);
        if (i0 == i1) throw new IllegalArgumentException("extremeIndex entries must be different");

        if (K == 2) {
            // only possible labeled move: swap which state is 0 vs 1
            extremeIndex.setParameterValue(0, i1);
            extremeIndex.setParameterValue(1, i0);

            // optional hard enforcement
            alphaRates.setParameterValue(i1, 0.0);
            alphaRates.setParameterValue(i0, 1.0);

        } else {

            // choose which labeled extreme to move: 0-label (slot 0) or 1-label (slot 1)
            final boolean moveOneLabel = MathUtils.nextBoolean();
            final int slot = moveOneLabel ? 1 : 0;
            final int i = moveOneLabel ? i1 : i0;
            final int other = moveOneLabel ? i0 : i1;

            // choose destination "k" uniformly among indices excluding {i, other} without rejection
            final int r = MathUtils.nextInt(K - 2);
            int k = r;
            if (k >= Math.min(i, other)) k++;
            if (k >= Math.max(i, other)) k++;

            // relabel: move chosen extreme label to k; k becomes old i
            extremeIndex.setParameterValue(slot, k);
            // keep the other slot unchanged (already equals other)

            double kValue = alphaRates.getParameterValue(k);
            double iValue = alphaRates.getParameterValue(i);

            alphaRates.setParameterValue(k, iValue);
            alphaRates.setParameterValue(i, kValue);

            if (tol >= 0.0) {
                final double a = alphaRates.getParameterValue(i);
                if (a < -tol || a > 1.0 + tol || Double.isNaN(a)) {
                    throw new IllegalStateException("alphaRates out of [0,1] at former extreme i=" + i + ": " + a);
                }
            }
        }

        alphaRates.fireParameterChangedEvent();
        extremeIndex.fireParameterChangedEvent();

        return 0.0;
    }

    private static int toIndexStrict(final double x, final String name) {
        final long r = Math.round(x);
        if (Math.abs(x - r) > 1e-9) {
            throw new IllegalArgumentException(name + " must be integer-like but found " + x);
        }
        if (r < Integer.MIN_VALUE || r > Integer.MAX_VALUE) {
            throw new IllegalArgumentException(name + " out of int range: " + x);
        }
        return (int) r;
    }
}