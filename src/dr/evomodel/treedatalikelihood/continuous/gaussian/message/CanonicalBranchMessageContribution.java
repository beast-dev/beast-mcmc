package dr.evomodel.treedatalikelihood.continuous.gaussian.message;

/**
 * Local canonical adjoints for a single branch factor
 * {@code p(x_{t+1} | x_t)} in block-canonical form.
 *
 * <p>For a branch factor written as
 * <pre>
 *   log p(z) = -0.5 z^T J z + h^T z - g,
 *   z = [x_t, x_{t+1}],
 * </pre>
 * Fisher's identity gives the local expected contribution under the branch posterior:
 * <pre>
 *   dL/dJ = -0.5 E[z z^T | Y]
 *   dL/dh =      E[z | Y]
 *   dL/dg =     -1
 * </pre>
 *
 * <p>Each block gradient matrix is stored row-major in a flat array of length
 * {@code dim * dim}: element {@code (i, j)} is at index {@code i * dim + j}.
 */
public final class CanonicalBranchMessageContribution {

    public final double[] dLogL_dPrecisionXX;
    public final double[] dLogL_dPrecisionXY;
    public final double[] dLogL_dPrecisionYX;
    public final double[] dLogL_dPrecisionYY;
    public final double[] dLogL_dInformationX;
    public final double[] dLogL_dInformationY;
    public double dLogL_dLogNormalizer;

    private final int dimension;

    public CanonicalBranchMessageContribution(final int dimension) {
        this.dimension = dimension;
        dLogL_dPrecisionXX = new double[dimension * dimension];
        dLogL_dPrecisionXY = new double[dimension * dimension];
        dLogL_dPrecisionYX = new double[dimension * dimension];
        dLogL_dPrecisionYY = new double[dimension * dimension];
        dLogL_dInformationX = new double[dimension];
        dLogL_dInformationY = new double[dimension];
        dLogL_dLogNormalizer = -1.0;
    }

    public int getDimension() {
        return dimension;
    }
}
