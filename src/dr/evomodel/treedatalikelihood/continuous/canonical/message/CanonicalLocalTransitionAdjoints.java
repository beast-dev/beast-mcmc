package dr.evomodel.treedatalikelihood.continuous.canonical.message;

/**
 * Local branch adjoints for a Gaussian transition in moment coordinates.
 *
 * <p>These are the transition-space sensitivities induced by one canonical
 * branch factor and its local posterior contribution:
 * {@code dL/dF}, {@code dL/df}, and {@code dL/dOmega}.
 *
 * <p>Matrix fields are stored row-major in flat arrays of length {@code dim * dim}:
 * element {@code (i, j)} is at index {@code i * dim + j}.
 */
public final class CanonicalLocalTransitionAdjoints {

    public final double[] dLogL_dF;
    public final double[] dLogL_df;
    public final double[] dLogL_dOmega;

    private final int dimension;

    public CanonicalLocalTransitionAdjoints(final int dimension) {
        this.dimension = dimension;
        this.dLogL_dF = new double[dimension * dimension];
        this.dLogL_df = new double[dimension];
        this.dLogL_dOmega = new double[dimension * dimension];
    }

    public int getDimension() {
        return dimension;
    }
}
