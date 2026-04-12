package dr.inference.timeseries.representation;

/**
 * Canonical-form representation of a Gaussian state density.
 *
 * <p>The density is represented as
 * <pre>
 *   p(x) = exp(-0.5 x^T J x + h^T x - g)
 * </pre>
 * where {@code J} is the precision matrix, {@code h} is the information vector,
 * and {@code g} stores the normalizing constant contribution.
 */
public final class CanonicalGaussianState {

    public final double[][] precision;
    public final double[] information;
    public double logNormalizer;

    public CanonicalGaussianState(final int dimension) {
        if (dimension < 1) {
            throw new IllegalArgumentException("dimension must be positive");
        }
        this.precision = new double[dimension][dimension];
        this.information = new double[dimension];
        this.logNormalizer = 0.0;
    }

    public int getDimension() {
        return information.length;
    }
}
