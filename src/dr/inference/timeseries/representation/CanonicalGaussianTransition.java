package dr.inference.timeseries.representation;

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
 */
public final class CanonicalGaussianTransition {

    public final double[][] precisionXX;
    public final double[][] precisionXY;
    public final double[][] precisionYX;
    public final double[][] precisionYY;
    public final double[] informationX;
    public final double[] informationY;
    public double logNormalizer;

    public CanonicalGaussianTransition(final int dimension) {
        if (dimension < 1) {
            throw new IllegalArgumentException("dimension must be positive");
        }
        this.precisionXX = new double[dimension][dimension];
        this.precisionXY = new double[dimension][dimension];
        this.precisionYX = new double[dimension][dimension];
        this.precisionYY = new double[dimension][dimension];
        this.informationX = new double[dimension];
        this.informationY = new double[dimension];
        this.logNormalizer = 0.0;
    }

    public int getDimension() {
        return informationX.length;
    }
}
