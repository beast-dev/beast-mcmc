//package dr.evomodel.treedatalikelihood.continuous.backprop;
//
//import org.ejml.data.DenseMatrix64F;
//
package dr.evomodel.treedatalikelihood.continuous.backprop;

import org.ejml.data.DenseMatrix64F;

import java.util.Arrays;

/**
 * Solves block Lyapunov systems of the form
 *
 *     D * Sigma + Sigma * D^T = V
 *
 * where D is encoded in compressed tridiagonal form
 *
 *     [diag(0..dim-1) | upper(0..dim-2) | lower(0..dim-2)]
 *
 * and the admissible block structure consists only of contiguous 1x1 and 2x2 blocks.
 *
 * This class is intentionally independent of any basis-change matrices.
 * It only knows about the raw entries of D in the block space.
 *
 * Not thread-safe: owns reusable workspaces.
 */
public final class BlockDiagonalLyapunovSolver {

    private static final double EPS = 1e-12;

    private final int dim;
    private final int upperOffset;
    private final int lowerOffset;

    /**
     * Cached structure inferred from the first parameter vector.
     * This is valid only if the 1x1/2x2 block pattern stays fixed.
     */
    private BlockStructure structure;

    // Reusable tiny dense workspaces.
    private final double[] a3 = new double[9];
    private final double[] b3 = new double[3];
    private final double[] x3 = new double[3];

    private final double[] a4 = new double[16];
    private final double[] b4 = new double[4];
    private final double[] x4 = new double[4];

    public BlockDiagonalLyapunovSolver(final int dim) {
        if (dim <= 0) {
            throw new IllegalArgumentException("dim must be positive");
        }
        this.dim = dim;
        this.upperOffset = dim;
        this.lowerOffset = dim + dim - 1;
    }

    public BlockDiagonalLyapunovSolver(final int dim,
                                       final int[] blockStarts,
                                       final int[] blockSizes) {
        this(dim);
        if (blockStarts == null || blockSizes == null) {
            throw new IllegalArgumentException("blockStarts and blockSizes must be non-null");
        }
        if (blockStarts.length != blockSizes.length) {
            throw new IllegalArgumentException("blockStarts and blockSizes must have the same length");
        }
        final int[] startsCopy = Arrays.copyOf(blockStarts, blockStarts.length);
        final int[] sizesCopy = Arrays.copyOf(blockSizes, blockSizes.length);
        validateStructure(startsCopy, sizesCopy);
        this.structure = new BlockStructure(startsCopy.length, startsCopy, sizesCopy);
    }

    public void solve(final double[] blockDParams,
                      final DenseMatrix64F rhs,
                      final DenseMatrix64F sigmaOut) {
        if (!trySolve(blockDParams, rhs, sigmaOut)) {
            // Parameters are in non-stationary region during leapfrog trajectory.
            // Signal -Infinity to the likelihood by filling with +Infinity.
            Arrays.fill(sigmaOut.data, Double.POSITIVE_INFINITY);
            // This will propagate to -Infinity log-likelihood, causing HMC to
            // reject the trajectory and reverse momentum correctly.
        }
    }

    /**
     * Non-throwing hot path for HMC/integration code.
     *
     * @return {@code true} when the Lyapunov system was solved successfully,
     * {@code false} when the current block parameters imply a singular or
     * non-stationary system.
     */
    public boolean trySolve(final double[] blockDParams,
                            final DenseMatrix64F rhs,
                            final DenseMatrix64F sigmaOut) {
        checkSquare(rhs, "rhs");
        checkSquare(sigmaOut, "sigmaOut");

        final BlockStructure blocks = ensureStructure(blockDParams);
        final double[] rhsData = rhs.data;
        final double[] outData = sigmaOut.data;

        Arrays.fill(outData, 0.0);

        for (int bi = 0; bi < blocks.count; bi++) {
            final int i0 = blocks.starts[bi];
            final int is = blocks.sizes[bi];
            for (int bj = 0; bj < blocks.count; bj++) {
                final int j0 = blocks.starts[bj];
                final int js = blocks.sizes[bj];
                if (!solvePair(blockDParams, rhsData, outData, i0, is, j0, js)) {
                    return false;
                }
            }
        }
        return true;
    }


    public void solving(final double[] blockDParams,
                      final DenseMatrix64F rhs,
                      final DenseMatrix64F sigmaOut) {
        if (!trySolve(blockDParams, rhs, sigmaOut)) {
            throw new LyapunovNonStationaryException("Singular or non-stationary block Lyapunov system");
        }
    }

    public void refreshBlockStructure(final double[] blockDParams) {
        this.structure = inferBlockStructure(blockDParams);
    }

    private BlockStructure ensureStructure(final double[] blockDParams) {
        BlockStructure s = structure;
        if (s == null) {
            s = inferBlockStructure(blockDParams);
            structure = s;
        }
        return s;
    }

    private BlockStructure inferBlockStructure(final double[] p) {
        final int[] startsTmp = new int[dim];
        final int[] sizesTmp = new int[dim];

        int count = 0;
        int i = 0;
        while (i < dim) {
            startsTmp[count] = i;
            sizesTmp[count] = isTwoByTwoBlock(p, i) ? 2 : 1;
            i += sizesTmp[count];
            count++;
        }

        final int[] starts = Arrays.copyOf(startsTmp, count);
        final int[] sizes = Arrays.copyOf(sizesTmp, count);
        return new BlockStructure(count, starts, sizes);
    }

    private void validateStructure(final int[] starts, final int[] sizes) {
        int next = 0;
        for (int i = 0; i < starts.length; ++i) {
            final int start = starts[i];
            final int size = sizes[i];
            if (start != next) {
                throw new IllegalArgumentException("Invalid block structure: expected block start " + next + " but found " + start);
            }
            if (size != 1 && size != 2) {
                throw new IllegalArgumentException("Only block sizes 1 and 2 are supported");
            }
            next += size;
        }
        if (next != dim) {
            throw new IllegalArgumentException("Block structure does not cover dimension " + dim);
        }
    }

    private boolean isTwoByTwoBlock(final double[] p, final int i) {
        if (i >= dim - 1) {
            return false;
        }
        return Math.abs(p[upperOffset + i]) >= EPS || Math.abs(p[lowerOffset + i]) >= EPS;
    }

    private boolean solvePair(final double[] p,
                              final double[] rhs,
                              final double[] out,
                              final int i0,
                              final int is,
                              final int j0,
                              final int js) {

        if (is == 1) {
            if (js == 1) {
                return solve11(p, rhs, out, i0, j0);
            } else {
                return solve12(p, rhs, out, i0, j0);
            }
        }

        if (js == 1) {
            return solve21(p, rhs, out, i0, j0);
        } else {
            return solve22(p, rhs, out, i0, j0);
        }
    }

    private boolean solve11(final double[] p,
                            final double[] rhs,
                            final double[] out,
                            final int i,
                            final int j) {

        final double denom = p[i] + p[j];
        if (!isNonSingular(denom)) {
            return false;
        }
        out[index(i, j)] = rhs[index(i, j)] / denom;
        return true;
    }

    private boolean solve12(final double[] p,
                            final double[] rhs,
                            final double[] out,
                            final int i,
                            final int j) {

        final double d = p[i];
        final double e00 = p[j];
        final double e01 = p[upperOffset + j];
        final double e10 = p[lowerOffset + j];
        final double e11 = p[j + 1];

        final int idx0 = index(i, j);
        final int idx1 = idx0 + 1;

        final double a11 = d + e00;
        final double a12 = e01;
        final double a21 = e10;
        final double a22 = d + e11;

        final double det = a11 * a22 - a12 * a21;
        if (!isNonSingular(det)) {
            return false;
        }

        final double v0 = rhs[idx0];
        final double v1 = rhs[idx1];

        out[idx0] = (v0 * a22 - v1 * a12) / det;
        out[idx1] = (a11 * v1 - a21 * v0) / det;
        return true;
    }

    private boolean solve21(final double[] p,
                            final double[] rhs,
                            final double[] out,
                            final int i,
                            final int j) {

        final double d00 = p[i];
        final double d01 = p[upperOffset + i];
        final double d10 = p[lowerOffset + i];
        final double d11 = p[i + 1];
        final double e = p[j];

        final int idx0 = index(i, j);
        final int idx1 = index(i + 1, j);

        final double a11 = d00 + e;
        final double a12 = d01;
        final double a21 = d10;
        final double a22 = d11 + e;

        final double det = a11 * a22 - a12 * a21;
        if (!isNonSingular(det)) {
            return false;
        }

        final double v0 = rhs[idx0];
        final double v1 = rhs[idx1];

        out[idx0] = (v0 * a22 - v1 * a12) / det;
        out[idx1] = (a11 * v1 - a21 * v0) / det;
        return true;
    }

    private boolean solve22(final double[] p,
                            final double[] rhs,
                            final double[] out,
                            final int i,
                            final int j) {

        final double d00 = p[i];
        final double d01 = p[upperOffset + i];
        final double d10 = p[lowerOffset + i];
        final double d11 = p[i + 1];

        final double e00 = p[j];
        final double e01 = p[upperOffset + j];
        final double e10 = p[lowerOffset + j];
        final double e11 = p[j + 1];

        final int ij00 = index(i, j);
        final int ij01 = ij00 + 1;
        final int ij10 = index(i + 1, j);
        final int ij11 = ij10 + 1;

        if (i == j) {
            // Symmetric 2x2 diagonal block.
            a3[0] = 2.0 * d00;  a3[1] = 2.0 * d01;  a3[2] = 0.0;
            a3[3] = d10;        a3[4] = d00 + d11; a3[5] = d01;
            a3[6] = 0.0;        a3[7] = 2.0 * d10; a3[8] = 2.0 * d11;

            b3[0] = rhs[ij00];
            b3[1] = rhs[ij01];
            b3[2] = rhs[ij11];

            if (!solveSmallSystem(3, a3, b3, x3)) {
                return false;
            }

            out[ij00] = x3[0];
            out[ij01] = x3[1];
            out[ij10] = x3[1];
            out[ij11] = x3[2];
            return true;
        }
        final boolean leftEqualDiagonal =
                Math.abs(d00 - d11) < EPS;

        final boolean rightEqualDiagonal =
                Math.abs(e00 - e11) < EPS;

        if (leftEqualDiagonal && rightEqualDiagonal) {

            final double[] h = b4;
            final double[] y = x4;

            h[0] = rhs[ij00];
            h[1] = rhs[ij01];
            h[2] = rhs[ij10];
            h[3] = rhs[ij11];

            if (!solveEqualDiagonal22(
                    d00, d01, d10,
                    e00, e01, e10,
                    h, y)) {
                return false;
            }

            out[ij00] = y[0];
            out[ij01] = y[1];
            out[ij10] = y[2];
            out[ij11] = y[3];

            return true;
        }

        a4[0]  = d00 + e00; a4[1]  = e01;       a4[2]  = d01;       a4[3]  = 0.0;
        a4[4]  = e10;       a4[5]  = d00 + e11; a4[6]  = 0.0;       a4[7]  = d01;
        a4[8]  = d10;       a4[9]  = 0.0;       a4[10] = d11 + e00; a4[11] = e01;
        a4[12] = 0.0;       a4[13] = d10;       a4[14] = e10;       a4[15] = d11 + e11;

        b4[0] = rhs[ij00];
        b4[1] = rhs[ij01];
        b4[2] = rhs[ij10];
        b4[3] = rhs[ij11];

        if (!solveSmallSystem(4, a4, b4, x4)) {
            return false;
        }

        out[ij00] = x4[0];
        out[ij01] = x4[1];
        out[ij10] = x4[2];
        out[ij11] = x4[3];
        return true;
    }

    /**
     * In-place Gaussian elimination with partial pivoting on tiny dense systems.
     * The coefficient array is overwritten.
     */
    private static boolean solveSmallSystem(final int n,
                                            final double[] a,
                                            final double[] b,
                                            final double[] xOut) {

        for (int col = 0; col < n; col++) {
            int pivot = col;
            double pivotAbs = Math.abs(a[col * n + col]);
            for (int row = col + 1; row < n; row++) {
                final double candidate = Math.abs(a[row * n + col]);
                if (candidate > pivotAbs) {
                    pivot = row;
                    pivotAbs = candidate;
                }
            }

            if (!isNonSingular(pivotAbs)) {
                return false;
            }
            if (pivot != col) {
                swapRows(a, n, col, pivot);
                swap(b, col, pivot);
            }

            final double diag = a[col * n + col];
            for (int row = col + 1; row < n; row++) {
                final double factor = a[row * n + col] / diag;
                if (factor == 0.0) {
                    continue;
                }
                a[row * n + col] = 0.0;
                for (int k = col + 1; k < n; k++) {
                    a[row * n + k] -= factor * a[col * n + k];
                }
                b[row] -= factor * b[col];
            }
        }

        for (int row = n - 1; row >= 0; row--) {
            double sum = b[row];
            for (int col = row + 1; col < n; col++) {
                sum -= a[row * n + col] * xOut[col];
            }
            xOut[row] = sum / a[row * n + row];
        }
        return true;
    }

    private static void swapRows(final double[] a, final int n, final int r1, final int r2) {
        final int base1 = r1 * n;
        final int base2 = r2 * n;
        for (int col = 0; col < n; col++) {
            final double tmp = a[base1 + col];
            a[base1 + col] = a[base2 + col];
            a[base2 + col] = tmp;
        }
    }

    private static void swap(final double[] x, final int i, final int j) {
        final double tmp = x[i];
        x[i] = x[j];
        x[j] = tmp;
    }

    private int index(final int row, final int col) {
        return row * dim + col;
    }

    private void checkSquare(final DenseMatrix64F m, final String name) {
        if (m.numRows != dim || m.numCols != dim) {
            throw new IllegalArgumentException(name + " must be " + dim + " x " + dim);
        }
    }

    private static boolean isNonSingular(final double value) {
        return Math.abs(value) >= EPS;
    }

    private static final class BlockStructure {
        final int count;
        final int[] starts;
        final int[] sizes;

        private BlockStructure(final int count, final int[] starts, final int[] sizes) {
            this.count = count;
            this.starts = starts;
            this.sizes = sizes;
        }
    }

    private static boolean solveEqualDiagonal22(final double ai,
                                                final double ui,
                                                final double vi,
                                                final double aj,
                                                final double uj,
                                                final double vj,
                                                final double[] h,
                                                final double[] y) {

        /*
         * Solve
         *
         *     D_i Y + Y D_j^T = H
         *
         * where
         *
         *     D_i = [ ai  ui ]
         *           [ vi  ai ]
         *
         *     D_j = [ aj  uj ]
         *           [ vj  aj ]
         *
         * using row-major ordering:
         *
         *     y = [ y00, y01, y10, y11 ].
         */

        final double s = ai + aj;

        final double h00 = h[0];
        final double h01 = h[1];
        final double h10 = h[2];
        final double h11 = h[3];

        /*
         * With
         *
         *     B = [ s   uj ]
         *         [ vj  s  ],
         *
         * the 4x4 system is
         *
         *     B yTop + ui yBottom = hTop
         *     vi yTop + B yBottom = hBottom.
         *
         * Therefore
         *
         *     (B^2 - ui vi I) yTop    = B hTop    - ui hBottom
         *     (B^2 - ui vi I) yBottom = B hBottom - vi hTop.
         */

        final double ri2 = ui * vi;
        final double rj2 = uj * vj;

        final double a = s * s + rj2 - ri2;
        final double b = 2.0 * s * uj;
        final double c = 2.0 * s * vj;

        final double det = a * a - b * c;

        if (!isNonSingular(det)) {
            return false;
        }

        final double invDet = 1.0 / det;

        final double q0 = s * h00 + uj * h01 - ui * h10;
        final double q1 = vj * h00 + s * h01 - ui * h11;

        final double r0 = s * h10 + uj * h11 - vi * h00;
        final double r1 = vj * h10 + s * h11 - vi * h01;

        y[0] = (a * q0 - b * q1) * invDet;
        y[1] = (a * q1 - c * q0) * invDet;
        y[2] = (a * r0 - b * r1) * invDet;
        y[3] = (a * r1 - c * r0) * invDet;

        return true;
    }

}