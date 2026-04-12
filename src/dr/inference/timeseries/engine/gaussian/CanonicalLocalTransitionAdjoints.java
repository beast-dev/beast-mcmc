package dr.inference.timeseries.engine.gaussian;

/**
 * Local branch adjoints for a Gaussian transition in moment coordinates.
 *
 * <p>These are the transition-space sensitivities induced by one canonical
 * branch factor and its local posterior contribution:
 * {@code dL/dF}, {@code dL/df}, and {@code dL/dOmega}.
 */
public final class CanonicalLocalTransitionAdjoints {

    public final double[][] dLogL_dF;
    public final double[] dLogL_df;
    public final double[][] dLogL_dOmega;

    public CanonicalLocalTransitionAdjoints(final int dimension) {
        this.dLogL_dF = new double[dimension][dimension];
        this.dLogL_df = new double[dimension];
        this.dLogL_dOmega = new double[dimension][dimension];
    }
}
