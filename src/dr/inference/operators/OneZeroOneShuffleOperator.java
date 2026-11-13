package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.math.MathUtils;

/**
 * Swaps one extreme reward label (0 or 1) with one non-extreme reward label
 * inside rewardRatesMapping.
 *
 * New parametrization:
 *
 *   rewardRatesValues[0] = 0
 *   rewardRatesValues[1] = 1
 *   rewardRatesValues[2..K-1] = internal continuous reward values
 *
 *   rewardRatesMapping[state] = reward-label index
 *
 * atomIndex indexes STATES, not reward labels, so it is not changed here.
 *
 * This operator changes which states are associated with the two extreme
 * reward labels 0 and 1.
 *
 * @author Filippo Monti
 */
public final class OneZeroOneShuffleOperator extends SimpleMCMCOperator {

    private final Parameter rewardRatesValues;    // length K
    private final Parameter rewardRatesMapping;   // length = number of states
    private final double tol;

    public OneZeroOneShuffleOperator(final Parameter rewardRatesValues,
                                     final Parameter rewardRatesMapping,
                                     final double weight,
                                     final double tol) {
        if (rewardRatesValues == null) {
            throw new IllegalArgumentException("rewardRatesValues must be non-null");
        }
        if (rewardRatesMapping == null) {
            throw new IllegalArgumentException("rewardRatesMapping must be non-null");
        }
        if (rewardRatesValues.getDimension() < 2) {
            throw new IllegalArgumentException("rewardRatesValues dimension must be >= 2");
        }
        if (tol < 0.0) {
            throw new IllegalArgumentException("tol must be >= 0");
        }

        this.rewardRatesValues = rewardRatesValues;
        this.rewardRatesMapping = rewardRatesMapping;
        this.tol = tol;

        validateRewardRatesValues();

        setWeight(weight);
    }

    @Override
    public String getOperatorName() {
        return "oneZeroOneShuffleOperator(" + rewardRatesMapping.getParameterName() + ")";
    }

    @Override
    public double doOperation() {

        final int K = rewardRatesValues.getDimension();

        localMove(K);

        rewardRatesMapping.fireParameterChangedEvent();

        return 0.0;
    }

    private void localMove(final int K) {

        if (K == 2) {
            // only labels 0 and 1 exist, so just swap them in the mapping
            remapRewardLabelsBySwap(0, 1);
            return;
        }

        // choose which extreme label to move: 0 or 1
        final boolean moveOneLabel = MathUtils.nextBoolean();
        final int i = moveOneLabel ? 1 : 0;

        // choose k uniformly from {2, ..., K-1}
        final int k = 2 + MathUtils.nextInt(K - 2);

        remapRewardLabelsBySwap(i, k);

        if (tol >= 0.0) {
            validateRewardRatesValues();
            validateRewardRatesMapping(K);
        }
    }

    /**
     * Swap reward labels i and k everywhere inside rewardRatesMapping.
     */
    private void remapRewardLabelsBySwap(final int i, final int k) {

        final int K = rewardRatesValues.getDimension();
        if (i < 0 || i >= K) {
            throw new IllegalArgumentException("swap index i out of range: " + i);
        }
        if (k < 0 || k >= K) {
            throw new IllegalArgumentException("swap index k out of range: " + k);
        }
        if (i == k) {
            return;
        }

        final int nStates = rewardRatesMapping.getDimension();
        for (int s = 0; s < nStates; s++) {
            final int r = toIndexStrict(
                    rewardRatesMapping.getParameterValue(s),
                    "rewardRatesMapping[" + s + "]"
            );

            if (r == i) {
                rewardRatesMapping.setParameterValueQuietly(s, k);
            } else if (r == k) {
                rewardRatesMapping.setParameterValueQuietly(s, i);
            }
        }
    }

    private void validateRewardRatesValues() {
        final int K = rewardRatesValues.getDimension();

        if (K < 2) {
            throw new IllegalArgumentException("rewardRatesValues must have dimension at least 2");
        }

        final double v0 = rewardRatesValues.getParameterValue(0);
        final double v1 = rewardRatesValues.getParameterValue(1);

        if (Math.abs(v0 - 0.0) > tol) {
            throw new IllegalArgumentException(
                    "rewardRatesValues[0] must be 0.0 but found " + v0
            );
        }
        if (Math.abs(v1 - 1.0) > tol) {
            throw new IllegalArgumentException(
                    "rewardRatesValues[1] must be 1.0 but found " + v1
            );
        }

        for (int j = 0; j < K; j++) {
            final double a = rewardRatesValues.getParameterValue(j);
            if (Double.isNaN(a) || a < -tol || a > 1.0 + tol) {
                throw new IllegalArgumentException(
                        "rewardRatesValues out of [0,1] at index " + j + ": " + a
                );
            }
        }
    }

    private void validateRewardRatesMapping(final int K) {
        final int nStates = rewardRatesMapping.getDimension();
        for (int s = 0; s < nStates; s++) {
            final int r = toIndexStrict(
                    rewardRatesMapping.getParameterValue(s),
                    "rewardRatesMapping[" + s + "]"
            );
            if (r < 0 || r >= K) {
                throw new IllegalStateException(
                        "rewardRatesMapping[" + s + "] out of range: " + r
                );
            }
        }
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