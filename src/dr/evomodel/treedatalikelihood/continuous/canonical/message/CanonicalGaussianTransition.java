package dr.evomodel.treedatalikelihood.continuous.canonical.message;

/**
 * Canonical-form representation of a Gaussian transition density.
 *
 * <p>For a branch transition from {@code x_t} to {@code x_{t+1}}, the conditional
 * density is represented as
 * <pre>
 *   p(x_{t+1} | x_t)
 *     = exp(-0.5 z^T J z + h^T z - g),
 *   z = [x_t, x_{t+1}]
 * </pre>
 * with block precision matrix
 * <pre>
 *   J = [[J_xx, J_xy],
 *        [J_yx, J_yy]]
 * </pre>
 * and block information vector {@code [h_x, h_y]}.
 *
 * <p>Each block is stored row-major in a flat array of length {@code dim * dim}:
 * element {@code (i, j)} is at index {@code i * dim + j}.
 */
public final class CanonicalGaussianTransition {

    public final double[] precisionXX;
    public final double[] precisionXY;
    public final double[] precisionYX;
    public final double[] precisionYY;
    public final double[] informationX;
    public final double[] informationY;
    public double logNormalizer;

    private final int dimension;

    public CanonicalGaussianTransition(final int dimension) {
        if (dimension < 1) {
            throw new IllegalArgumentException("dimension must be positive");
        }
        this.dimension = dimension;
        this.precisionXX = new double[dimension * dimension];
        this.precisionXY = new double[dimension * dimension];
        this.precisionYX = new double[dimension * dimension];
        this.precisionYY = new double[dimension * dimension];
        this.informationX = new double[dimension];
        this.informationY = new double[dimension];
        this.logNormalizer = 0.0;
    }

    public int getDimension() {
        return dimension;
    }

    public void copyFrom(final CanonicalGaussianTransition source) {
        if (source.dimension != dimension) {
            throw new IllegalArgumentException(
                    "Transition dimension mismatch: " + source.dimension + " vs " + dimension);
        }
        final int d2 = dimension * dimension;
        System.arraycopy(source.informationX, 0, informationX, 0, dimension);
        System.arraycopy(source.informationY, 0, informationY, 0, dimension);
        System.arraycopy(source.precisionXX, 0, precisionXX, 0, d2);
        System.arraycopy(source.precisionXY, 0, precisionXY, 0, d2);
        System.arraycopy(source.precisionYX, 0, precisionYX, 0, d2);
        System.arraycopy(source.precisionYY, 0, precisionYY, 0, d2);
        logNormalizer = source.logNormalizer;
    }
}
