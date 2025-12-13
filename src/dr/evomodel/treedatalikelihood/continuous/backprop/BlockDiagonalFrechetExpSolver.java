package dr.evomodel.treedatalikelihood.continuous.backprop;

import org.ejml.data.DenseMatrix64F;
import java.util.Arrays;

/**
 * Block-diagonal Fréchet exponential solver in the D basis.
 *
 * Solves, for block-diagonal D,
 *
 *        D_i X_{ij} + X_{ij} D_j = E_{ij},
 *
 * for every pair of blocks (i,j).  This is Option A: the Fréchet equation
 * has identical structure to the Sylvester equation used by the Lyapunov
 * solver, except that D_j is *not transposed*.
 *
 * We reuse the exact block structure, 1×1, 1×2, 2×1, 2×2 logic, and
 * the same 3×3 and 4×4 internal solvers. Complexity O(k^2).
 *
 * This is the exponential Fréchet derivative *in the D basis*. The caller
 * is expected to wrap this with similarity transforms:
 *
 *      L_exp(S; E) = R · L_exp(D; R^{-1} E R ) · R^{-1}.
 *
 * Block structure is computed once and cached.
 */
public final class BlockDiagonalFrechetExpSolver {

    private static final double EPS = 1e-12;

    private final int dim;

    /** Cached block structure (1×1 or 2×2 blocks). */
    private BlockStructureInfo blockInfo = null;

    /** Reusable temporary arrays for 3×3 and 4×4 solvers. */
    private final double[] A3 = new double[9];
    private final double[] b3 = new double[3];
    private final double[] x3 = new double[3];

    private final double[] A4 = new double[16];
    private final double[] b4 = new double[4];
    private final double[] x4 = new double[4];

    public BlockDiagonalFrechetExpSolver(int dim) {
        this.dim = dim;
    }

    // ---------------------------------------------------------------------
    // Public solve
    // ---------------------------------------------------------------------

    /**
     * Solve, for block-diagonal D,
     *
     *     D X + X D = E
     *
     * in the D basis.
     *
     * blockDParams is in compressed format:
     *   [ d_0..d_{d-1}, u_0..u_{d-2}, l_0..l_{d-2} ]
     */
    public void solve(double[] blockDParams,
                      DenseMatrix64F E,
                      DenseMatrix64F XOut)
    {
        if (blockInfo == null)
            blockInfo = identifyBlockStructure(blockDParams, dim);

        Arrays.fill(XOut.data, 0.0);

        final double[] eData = E.data;
        final double[] xData = XOut.data;

        for (int bi = 0; bi < blockInfo.numBlocks; bi++) {
            int iStart = blockInfo.offsets[bi];
            int iSize  = blockInfo.sizes[bi];

            for (int bj = 0; bj < blockInfo.numBlocks; bj++) {
                int jStart = blockInfo.offsets[bj];
                int jSize  = blockInfo.sizes[bj];

                solveBlockFréchet(
                        blockDParams,
                        eData, xData,
                        iStart, iSize,
                        jStart, jSize
                );
            }
        }
    }


    // ---------------------------------------------------------------------
    // Block structure
    // ---------------------------------------------------------------------

    private static final class BlockStructureInfo {
        final int numBlocks;
        final int[] offsets;
        final int[] sizes;
        BlockStructureInfo(int n, int[] off, int[] sz) {
            numBlocks = n;
            offsets   = off;
            sizes     = sz;
        }
    }

    /**
     * Identical logic to the Lyapunov solver:
     * A 2×2 block occurs when BOTH u_i and l_i are nonzero.
     */
    private static BlockStructureInfo identifyBlockStructure(double[] blockDParams, int dim) {
        int upperOffset = dim;
        int lowerOffset = dim + (dim - 1);

        int[] tmpOff = new int[dim];
        int[] tmpSz  = new int[dim];
        int count = 0;

        int i = 0;
        while (i < dim) {
            tmpOff[count] = i;
            if (i == dim - 1) {
                tmpSz[count] = 1;
                i++;
            } else {
                double u = blockDParams[upperOffset + i];
                double l = blockDParams[lowerOffset + i];
                if (Math.abs(u) < EPS && Math.abs(l) < EPS) {
                    tmpSz[count] = 1;
                    i++;
                } else {
                    tmpSz[count] = 2;
                    i += 2;
                }
            }
            count++;
        }

        int[] off = new int[count];
        int[] sz  = new int[count];
        System.arraycopy(tmpOff, 0, off, 0, count);
        System.arraycopy(tmpSz,  0, sz,  0, count);

        return new BlockStructureInfo(count, off, sz);
    }

    // ---------------------------------------------------------------------
    // Block dispatcher
    // ---------------------------------------------------------------------

    private void solveBlockFréchet(double[] blockDParams,
                                   double[] eData,
                                   double[] xData,
                                   int iStart, int iSize,
                                   int jStart, int jSize)
    {
        if (iSize == 1 && jSize == 1) {
            solve1x1(blockDParams, eData, xData, iStart, jStart);
        }
        else if (iSize == 1 && jSize == 2) {
            solve1x2(blockDParams, eData, xData, iStart, jStart);
        }
        else if (iSize == 2 && jSize == 1) {
            solve2x1(blockDParams, eData, xData, iStart, jStart);
        }
        else { // 2×2 × 2×2
            solve2x2(blockDParams, eData, xData, iStart, jStart);
        }
    }

    // ---------------------------------------------------------------------
    // 1×1 × 1×1
    // ---------------------------------------------------------------------

    /** Solve (d_i + d_j) x = e. */
    private void solve1x1(double[] p, double[] e, double[] x, int i, int j) {
        double di = p[i];
        double dj = p[j];

        double denom = di + dj;
        if (Math.abs(denom) < EPS)
            throw new RuntimeException("Singular 1×1 Fréchet block at (" + i + "," + j + ")");

        int idx = i * dim + j;
        x[idx] = e[idx] / denom;
    }

    // ---------------------------------------------------------------------
    // 1×1 × 2×2
    // ---------------------------------------------------------------------

    private void solve1x2(double[] p, double[] e, double[] x, int i, int j) {
        int uOff = dim;
        int lOff = dim + (dim - 1);

        double d  = p[i];
        double d0 = p[j];
        double d1 = p[uOff + j];
        double d2 = p[lOff + j];
        double d3 = p[j + 1];

        int idx0 = i * dim + j;
        int idx1 = idx0 + 1;

        double e0 = e[idx0];
        double e1 = e[idx1];

        double a11 = d + d0;
        double a12 = d2;
        double a21 = d1;
        double a22 = d + d3;

        double det = a11 * a22 - a12 * a21;
        if (Math.abs(det) < EPS)
            throw new RuntimeException("Singular 1×2 block at (" + i + "," + j + ")");

        x[idx0] = (e0 * a22 - e1 * a12) / det;
        x[idx1] = (a11 * e1 - a21 * e0) / det;
    }

    // ---------------------------------------------------------------------
    // 2×2 × 1×1
    // ---------------------------------------------------------------------

    private void solve2x1(double[] p, double[] e, double[] x, int i, int j) {
        int uOff = dim;
        int lOff = dim + (dim - 1);

        double d0 = p[i];
        double d1 = p[uOff + i];
        double d2 = p[lOff + i];
        double d3 = p[i + 1];
        double d  = p[j];

        int idx0 = i * dim + j;
        int idx1 = (i + 1) * dim + j;

        double e0 = e[idx0];
        double e1 = e[idx1];

        double a11 = d0 + d;
        double a12 = d1;
        double a21 = d2;
        double a22 = d3 + d;

        double det = a11 * a22 - a12 * a21;
        if (Math.abs(det) < EPS)
            throw new RuntimeException("Singular 2×1 block at (" + i + "," + j + ")");

        x[idx0] = (e0 * a22 - e1 * a12) / det;
        x[idx1] = (a11 * e1 - a21 * e0) / det;
    }

    // ---------------------------------------------------------------------
    // 2×2 × 2×2
    // ---------------------------------------------------------------------

    private void solve2x2(double[] p, double[] e, double[] x, int i, int j)
    {
        int uOff = dim;
        int lOff = dim + (dim - 1);

        // Di block
        double d0 = p[i];
        double d1 = p[uOff + i];
        double d2 = p[lOff + i];
        double d3 = p[i + 1];

        // Dj block
        double e0b = p[j];
        double e1b = p[uOff + j];
        double e2b = p[lOff + j];
        double e3b = p[j + 1];

        // Indices in full matrices
        int idx00 = i       * dim + j;
        int idx01 = idx00 + 1;
        int idx10 = (i + 1) * dim + j;
        int idx11 = idx10 + 1;

        double E00 = e[idx00];
        double E01 = e[idx01];
        double E10 = e[idx10];
        double E11 = e[idx11];

        // diagonal block (i == j) → symmetric 3×3
        if (i == j) {
            double[] A = A3;
            double[] B = b3;
            double[] X = x3;

            // Equations for variables (x00, x01, x11)
            A[0] = 2*d0;   A[1] = 2*d1;   A[2] = 0.0;
            A[3] = d2;     A[4] = d0+d3; A[5] = d1;
            A[6] = 0.0;    A[7] = 2*d2;   A[8] = 2*d3;

            B[0] = E00;
            B[1] = E01;
            B[2] = E11;

            solve3x3(A, B, X);

            double x00 = X[0];
            double x01 = X[1];
            double x11 = X[2];

            x[idx00] = x00;
            x[idx01] = x01;
            x[idx10] = x01;
            x[idx11] = x11;
        }
        else {
            double[] A = A4;
            double[] B = b4;
            double[] X = x4;

            // Fill 4×4 system (no transpose on Dj!)
            A[0]  = d0 + e0b;   A[1]  = d1;        A[2]  = e2b;       A[3]  = 0.0;
            A[4]  = d2;         A[5]  = d0 + e3b;  A[6]  = 0.0;       A[7]  = e2b;
            A[8]  = e1b;        A[9]  = 0.0;       A[10] = d3 + e0b;  A[11] = d1;
            A[12] = 0.0;        A[13] = e1b;       A[14] = d2;        A[15] = d3 + e3b;

            B[0] = E00;
            B[1] = E01;
            B[2] = E10;
            B[3] = E11;

            solve4x4(A, B, X);

            x[idx00] = X[0];
            x[idx01] = X[1];
            x[idx10] = X[2];
            x[idx11] = X[3];
        }
    }

    // ---------------------------------------------------------------------
    // Small system solvers (identical to Lyapunov solver)
    // ---------------------------------------------------------------------

    private void solve3x3(double[] A, double[] b, double[] x) {
        double a00=A[0], a01=A[1], a02=A[2];
        double a10=A[3], a11=A[4], a12=A[5];
        double a20=A[6], a21=A[7], a22=A[8];

        double b0=b[0], b1=b[1], b2=b[2];

        double det =
                a00*(a11*a22 - a12*a21)
                        - a01*(a10*a22 - a12*a20)
                        + a02*(a10*a21 - a11*a20);

        if (Math.abs(det) < EPS)
            throw new RuntimeException("Singular 3×3 Fréchet system");

        double det0 =
                b0*(a11*a22 - a12*a21)
                        - a01*(b1*a22 - a12*b2)
                        + a02*(b1*a21 - a11*b2);

        double det1 =
                a00*(b1*a22 - a12*b2)
                        - b0*(a10*a22 - a12*a20)
                        + a02*(a10*b2 - b1*a20);

        double det2 =
                a00*(a11*b2 - b1*a21)
                        - a01*(a10*b2 - b1*a20)
                        + b0*(a10*a21 - a11*a20);

        double inv = 1.0/det;

        x[0] = det0*inv;
        x[1] = det1*inv;
        x[2] = det2*inv;
    }

    private void solve4x4(double[] A, double[] b, double[] x) {

        // Forward elimination
        for (int col = 0; col < 4; col++) {

            // pivot
            int piv = col;
            double max = Math.abs(A[col*4+col]);
            for (int r = col+1; r < 4; r++) {
                double v = Math.abs(A[r*4+col]);
                if (v > max) { max = v; piv = r; }
            }

            if (max < EPS)
                throw new RuntimeException("Singular 4×4 Fréchet system");

            // row swap
            if (piv != col) {
                int c0 = col*4, cp = piv*4;
                for (int k = 0; k < 4; k++) {
                    double tmp = A[c0+k];
                    A[c0+k] = A[cp+k];
                    A[cp+k] = tmp;
                }
                double tb = b[col];
                b[col]    = b[piv];
                b[piv]    = tb;
            }

            // eliminate
            double pivot = A[col*4 + col];
            for (int r = col+1; r < 4; r++) {
                int rr = r*4;
                double fac = A[rr+col] / pivot;
                if (fac != 0.0) {
                    for (int k = col; k < 4; k++)
                        A[rr+k] -= fac * A[col*4+k];
                    b[r] -= fac * b[col];
                }
            }
        }

        // back substitution
        for (int r = 3; r >= 0; r--) {
            double sum = b[r];
            int rr = r*4;
            for (int c = r+1; c < 4; c++)
                sum -= A[rr+c] * x[c];
            double diag = A[rr+r];
            if (Math.abs(diag) < EPS)
                throw new RuntimeException("Singular 4×4 Fréchet system");
            x[r] = sum / diag;
        }
    }
}

