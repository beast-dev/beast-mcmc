package dr.evomodel.treedatalikelihood.continuous.gaussian;

/**
 * Canonical-form representation of a Gaussian state density.
 *
 * <p>The density is represented as
 * <pre>
 *   p(x) = exp(-0.5 x^T J x + h^T x - g)
 * </pre>
 * where {@code J} is the precision matrix, {@code h} is the information vector,
 * and {@code g} stores the normalizing constant contribution.
 *
 * <p>The precision matrix is stored row-major in a flat array of length
 * {@code dim * dim}: element {@code (i, j)} is at index {@code i * dim + j}.
 */
public final class CanonicalGaussianState {

    public final double[] precision;
    public final double[] information;
    public double logNormalizer;

    private final int dimension;

    public CanonicalGaussianState(final int dimension) {
        if (dimension < 1) {
            throw new IllegalArgumentException("dimension must be positive");
        }
        this.dimension = dimension;
        this.precision = new double[dimension * dimension];
        this.information = new double[dimension];
        this.logNormalizer = 0.0;
    }

    public int getDimension() {
        return dimension;
    }
}
