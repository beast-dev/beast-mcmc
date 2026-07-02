package dr.evomodel.branchratemodel;

import dr.inference.hmc.EmbeddedOrdinalParameter;
import dr.inference.model.Parameter;

/**
 * Decodes a single embedded categorical branch-state parameter.
 *
 * Category 0 is the continuous reward state. Categories 1..K are atomic
 * no-jump reward states 0..K-1.
 */
public final class RewardMixtureCategoryDecoder {

    public static final int CONTINUOUS_CATEGORY = 0;

    private final Parameter categoryParameter;
    private final Parameter cutParameter;
    private final int atomicStateCount;
    private EmbeddedOrdinalParameter embedding;

    public RewardMixtureCategoryDecoder(final Parameter categoryParameter,
                                        final Parameter cutParameter,
                                        final int atomicStateCount,
                                        final int expectedDimension) {
        if (categoryParameter == null) {
            throw new IllegalArgumentException("categoryParameter must be non-null");
        }
        if (cutParameter == null) {
            throw new IllegalArgumentException("cutParameter must be non-null");
        }
        if (atomicStateCount < 1) {
            throw new IllegalArgumentException("atomicStateCount must be positive");
        }
        if (categoryParameter.getDimension() != expectedDimension) {
            throw new IllegalArgumentException(
                    "categoryParameter dimension must be " + expectedDimension +
                            " but is " + categoryParameter.getDimension());
        }
        if (cutParameter.getDimension() != atomicStateCount + 2) {
            throw new IllegalArgumentException(
                    "cutParameter dimension must be atomicStateCount + 2 (" +
                            (atomicStateCount + 2) + ") but is " + cutParameter.getDimension());
        }

        this.categoryParameter = categoryParameter;
        this.cutParameter = cutParameter;
        this.atomicStateCount = atomicStateCount;
        refreshEmbedding();
    }

    public Parameter getCategoryParameter() {
        return categoryParameter;
    }

    public Parameter getCutParameter() {
        return cutParameter;
    }

    public int getCategoryCount() {
        return atomicStateCount + 1;
    }

    public void refreshEmbedding() {
        embedding = new EmbeddedOrdinalParameter(cutParameter.getParameterValues());
        if (embedding.getStateCount() != getCategoryCount()) {
            throw new IllegalArgumentException(
                    "Embedding state count must be " + getCategoryCount() +
                            " but is " + embedding.getStateCount());
        }
    }

    public int getCategoryForParameterIndex(final int parameterIndex) {
        if (parameterIndex < 0 || parameterIndex >= categoryParameter.getDimension()) {
            throw new IllegalArgumentException("Category parameter index out of range: " + parameterIndex);
        }
        try {
            return getCategoryForValue(categoryParameter.getParameterValue(parameterIndex));
        } catch (IllegalArgumentException e) {
            throw new IllegalArgumentException(
                    "Invalid embedded category value at parameter index " + parameterIndex +
                            ": " + categoryParameter.getParameterValue(parameterIndex), e);
        }
    }

    public int getCategoryForValue(final double value) {
        final int category = embedding.getStateIndex(value);
        if (category < 0 || category >= getCategoryCount()) {
            throw new IllegalArgumentException("Decoded category out of range: " + category);
        }
        return category;
    }

    public double getLowerCut(final int category) {
        return embedding.getLowerCut(category);
    }

    public double getUpperCut(final int category) {
        return embedding.getUpperCut(category);
    }

    public double getNextBoundary(final double value, final double direction) {
        return embedding.getNextBoundary(value, direction);
    }

    public boolean isAtomic(final int parameterIndex) {
        return isAtomicCategory(getCategoryForParameterIndex(parameterIndex));
    }

    public int getAtomicState(final int parameterIndex) {
        return getAtomicStateForCategory(getCategoryForParameterIndex(parameterIndex));
    }

    public static boolean isContinuousCategory(final int category) {
        return category == CONTINUOUS_CATEGORY;
    }

    public static boolean isAtomicCategory(final int category) {
        return category > CONTINUOUS_CATEGORY;
    }

    public static int getAtomicStateForCategory(final int category) {
        if (!isAtomicCategory(category)) {
            throw new IllegalArgumentException("Category " + category + " is not atomic");
        }
        return category - 1;
    }
}
