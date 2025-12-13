package dr.evomodel.treedatalikelihood.continuous.backprop;

import org.ejml.data.DenseMatrix64F;

import java.util.Arrays;

/**
 * Helper for solving block-diagonal Lyapunov equations of the form
 *
 *      D Σ + Σ Dᵀ = V
 *
 * where D is given in compressed block-tridiagonal form with 1×1 and 2×2
 * real blocks:
 *
 *   blockDParams = [ d_0, ..., d_{d-1},
 *                    u_0, ..., u_{d-2},
 *                    l_0, ..., l_{d-2} ]
 *
 * D has
 *   D_{ii}     = d_i,
 *   D_{i,i+1}  = u_i,
 *   D_{i+1,i}  = l_i,
 * and all other entries zero.  Contiguous non-zero pairs (u_i, l_i)
 * represent 2×2 blocks; all other positions are 1×1 blocks.
 *
 * This class:
 *   - Identifies the block structure of D once and caches it.
 *   - Solves the Lyapunov equation blockwise by iterating over all block
 *     pairs (i, j).
 *   - Each block–block equation is solved analytically using fixed-size
 *     1×1, 2×2, 3×3, or 4×4 linear systems, with O(1) cost per block pair.
 *   - With O(n) blocks in an n-dimensional system, the overall complexity
 *     is O(n²) in the process dimension. (counting the pairwise block solves)
 *
 * The block structure (i.e. the pattern of 1×1 versus 2×2 blocks) is assumed
 * to remain fixed for the lifetime of this solver instance, although the
 * numerical values in blockDParams may change during MCMC.
 */

/**
 * @author Filippo Monti
 */
public final class BlockDiagonalLyapunovSolver {

    private static final double EPS_BLOCK = 1e-12;

    private final int dim;

    // Cached block structure (computed on first use)
    private BlockStructureInfo blockStructureInfo = null;

    // Small reusable workspaces (avoid allocations in hot paths).
    // 3×3 system: A(3×3), b(3), x(3) in flat row-major layout.
    private final double[] lyapA3 = new double[9];
    private final double[] lyapB3 = new double[3];
    private final double[] lyapX3 = new double[3];

    // 4×4 system: A(4×4), b(4), x(4) in flat row-major layout.
    private final double[] lyapA4 = new double[16];
    private final double[] lyapB4 = new double[4];
    private final double[] lyapX4 = new double[4];

    public BlockDiagonalLyapunovSolver(int dim) {
        this.dim = dim;
    }

    /**
     * Solve the continuous Lyapunov equation
     *
     *     D Σ + Σ Dᵀ = V
     *
     * using the current block structure of D implied by blockDParams.
     *
     * @param blockDParams compressed block-tridiagonal parameters of D
     * @param V            right-hand side matrix (dim × dim)
     * @param Sigma        output matrix: solution Σ (dim × dim)
     */
    public void solve(double[] blockDParams,
                      DenseMatrix64F V,
                      DenseMatrix64F Sigma) {

        // Compute and cache block structure on first use.
        BlockStructureInfo blocks = blockStructureInfo;
        if (blocks == null) {
            blocks = identifyBlockStructure(blockDParams, dim);
            blockStructureInfo = blocks;
        }

        double[] vData = V.data;
        double[] sData = Sigma.data;

        // Zero-out Sigma
        Arrays.fill(sData, 0.0);

        // For each block pair (i,j), solve:
        //   Dᵢ Σᵢⱼ + Σᵢⱼ Dⱼᵀ = Vᵢⱼ
        for (int bi = 0; bi < blocks.numBlocks; bi++) {
            int iStart = blocks.offsets[bi];
            int iSize  = blocks.sizes[bi];

            for (int bj = 0; bj < blocks.numBlocks; bj++) {
                int jStart = blocks.offsets[bj];
                int jSize  = blocks.sizes[bj];

                solveBlockLyapunov(blockDParams,
                        vData, sData,
                        iStart, iSize,
                        jStart, jSize);
            }
        }
    }

    // ----------------------------------------------------------------------
    // Block structure
    // ----------------------------------------------------------------------

    private static final class BlockStructureInfo {
        final int numBlocks;
        final int[] offsets;  // starting index of each block
        final int[] sizes;    // size of each block (1 or 2)

        BlockStructureInfo(int numBlocks, int[] offsets, int[] sizes) {
            this.numBlocks = numBlocks;
            this.offsets = offsets;
            this.sizes = sizes;
        }
    }

    /**
     * Identify 1×1 and 2×2 blocks directly from blockDParams.
     *
     * blockDParams layout: [d₀...d_{n-1}, u₀...u_{n-2}, l₀...l_{n-2}]
     * where D[i,i]=d_i, D[i,i+1]=u_i, D[i+1,i]=l_i.
     * Contiguous non-zero (u_i, l_i) pairs define 2×2 blocks.
     */
    private static BlockStructureInfo identifyBlockStructure(double[] blockDParams, int dim) {
        int upperOffset = dim;
        int lowerOffset = dim + (dim - 1);

        int[] tempOffsets = new int[dim];
        int[] tempSizes   = new int[dim];
        int blockCount    = 0;

        int i = 0;
        while (i < dim) {
            tempOffsets[blockCount] = i;

            if (i == dim - 1) {
                // Last element is always 1×1
                tempSizes[blockCount] = 1;
                i++;
            } else {
                double u = blockDParams[upperOffset + i];
                double l = blockDParams[lowerOffset + i];

                if (Math.abs(u) < EPS_BLOCK && Math.abs(l) < EPS_BLOCK) {
                    // 1×1 block
                    tempSizes[blockCount] = 1;
                    i++;
                } else {
                    // 2×2 block
                    tempSizes[blockCount] = 2;
                    i += 2;
                }
            }
            blockCount++;
        }

        int[] offsets = new int[blockCount];
        int[] sizes   = new int[blockCount];
        System.arraycopy(tempOffsets, 0, offsets, 0, blockCount);
        System.arraycopy(tempSizes,   0, sizes,   0, blockCount);

        return new BlockStructureInfo(blockCount, offsets, sizes);
    }

    // ----------------------------------------------------------------------
    // Block-level Lyapunov solves
    // ----------------------------------------------------------------------

    /**
     * Solve a single block equation:
     *
     *     Dᵢ Σᵢⱼ + Σᵢⱼ Dⱼᵀ = Vᵢⱼ
     *
     * iStart, jStart are row/col indices; iSize, jSize ∈ {1,2}.
     *
     * Cases:
     *   1×1 × 1×1: σ = v / (dᵢ + dⱼ)
     *   1×1 × 2×2: 2×2 system
     *   2×2 × 1×1: 2×2 system
     *   2×2 × 2×2: 3×3 (diagonal) or 4×4 (off-diagonal) system
     */
    private void solveBlockLyapunov(double[] blockDParams,
                                    double[] vData,
                                    double[] sData,
                                    int iStart, int iSize,
                                    int jStart, int jSize) {

        if (iSize == 1 && jSize == 1) {
            // 1×1 × 1×1 case
            double di = blockDParams[iStart];
            double dj = blockDParams[jStart];

            int idx = iStart * dim + jStart;
            double v = vData[idx];

            double denom = di + dj;
            if (Math.abs(denom) < EPS_BLOCK) {
                throw new RuntimeException("Singular 1×1 Lyapunov block at (" + iStart + "," + jStart + ")");
            }

            sData[idx] = v / denom;

        } else if (iSize == 1 && jSize == 2) {
            solve1x2BlockLyapunov(blockDParams, vData, sData, iStart, jStart);

        } else if (iSize == 2 && jSize == 1) {
            solve2x1BlockLyapunov(blockDParams, vData, sData, iStart, jStart);

        } else if (iSize == 2 && jSize == 2) {
            solve2x2BlockLyapunov(blockDParams, vData, sData, iStart, jStart);
        }
    }

    /**
     * 1×1 × 2×2 block:
     *
     *   d * [σ₀, σ₁] + [σ₀, σ₁] * [d₀  d₁] = [v₀, v₁]
     *                                  [d₂  d₃]
     *
     * → 2×2 system:
     *   (d + d₀) σ₀ + d₂ σ₁ = v₀
     *   d₁ σ₀ + (d + d₃) σ₁ = v₁
     */
    private void solve1x2BlockLyapunov(double[] blockDParams,
                                       double[] vData,
                                       double[] sData,
                                       int i, int j) {

        final int upperOffset = dim;
        final int lowerOffset = dim + (dim - 1);

        double d  = blockDParams[i];
        double d0 = blockDParams[j];
        double d1 = blockDParams[upperOffset + j];
        double d2 = blockDParams[lowerOffset + j];
        double d3 = blockDParams[j + 1];

        int idx0 = i * dim + j;
        int idx1 = idx0 + 1;

        double v0 = vData[idx0];
        double v1 = vData[idx1];

        double a11 = d + d0;
        double a12 = d2;
        double a21 = d1;
        double a22 = d + d3;

        double det = a11 * a22 - a12 * a21;
        if (Math.abs(det) < EPS_BLOCK) {
            throw new RuntimeException("Singular Lyapunov block at (" + i + "," + j + ")");
        }

        double s0 = (v0 * a22 - v1 * a12) / det;
        double s1 = (a11 * v1 - a21 * v0) / det;

        sData[idx0] = s0;
        sData[idx1] = s1;
    }

    /**
     * 2×2 × 1×1 block:
     *
     *   [d₀  d₁] [σ₀] + [σ₀] d = [v₀]
     *   [d₂  d₃] [σ₁]   [σ₁]     [v₁]
     *
     * → 2×2 system:
     *   (d₀ + d) σ₀ + d₁ σ₁ = v₀
     *   d₂ σ₀ + (d₃ + d) σ₁ = v₁
     */
    private void solve2x1BlockLyapunov(double[] blockDParams,
                                       double[] vData,
                                       double[] sData,
                                       int i, int j) {

        final int upperOffset = dim;
        final int lowerOffset = dim + (dim - 1);

        double d0 = blockDParams[i];
        double d1 = blockDParams[upperOffset + i];
        double d2 = blockDParams[lowerOffset + i];
        double d3 = blockDParams[i + 1];
        double d  = blockDParams[j];

        int idx0 = i * dim + j;
        int idx1 = (i + 1) * dim + j;

        double v0 = vData[idx0];
        double v1 = vData[idx1];

        double a11 = d0 + d;
        double a12 = d1;
        double a21 = d2;
        double a22 = d3 + d;

        double det = a11 * a22 - a12 * a21;
        if (Math.abs(det) < EPS_BLOCK) {
            throw new RuntimeException("Singular Lyapunov block at (" + i + "," + j + ")");
        }

        double s0 = (v0 * a22 - v1 * a12) / det;
        double s1 = (a11 * v1 - a21 * v0) / det;

        sData[idx0] = s0;
        sData[idx1] = s1;
    }

    /**
     * 2×2 × 2×2 block.
     *
     * Diagonal case (i == j): Σ is symmetric, reduced 3×3 system.
     * Off-diagonal: full 4×4 system.
     */
    private void solve2x2BlockLyapunov(double[] blockDParams,
                                       double[] vData,
                                       double[] sData,
                                       int i, int j) {

        final int upperOffset = dim;
        final int lowerOffset = dim + (dim - 1);

        // Di block
        double d0 = blockDParams[i];
        double d1 = blockDParams[upperOffset + i];
        double d2 = blockDParams[lowerOffset + i];
        double d3 = blockDParams[i + 1];

        // Dj block
        double e0 = blockDParams[j];
        double e1 = blockDParams[upperOffset + j];
        double e2 = blockDParams[lowerOffset + j];
        double e3 = blockDParams[j + 1];

        int idx00 = i       * dim + j;
        int idx01 = idx00 + 1;
        int idx10 = (i + 1) * dim + j;
        int idx11 = idx10 + 1;

        double v0 = vData[idx00];
        double v1 = vData[idx01];
        double v2 = vData[idx10];
        double v3 = vData[idx11];

        if (i == j) {
            // Symmetric diagonal block: 3×3 system
            double[] A = lyapA3;
            double[] b = lyapB3;
            double[] x = lyapX3;

            // 3×3 system for [σ₀, σ₁, σ₃]ᵀ
            // A is row-major 3×3.
            A[0] = 2.0 * d0;  A[1] = 2.0 * d1;  A[2] = 0.0;
            A[3] = d2;        A[4] = d0 + d3;  A[5] = d1;
            A[6] = 0.0;       A[7] = 2.0 * d2; A[8] = 2.0 * d3;

            b[0] = v0;
            b[1] = v1;
            b[2] = v3;

            solve3x3InPlace(A, b, x);

            double s0 = x[0];
            double s1 = x[1];
            double s3 = x[2];

            sData[idx00] = s0;
            sData[idx01] = s1;
            sData[idx10] = s1;   // symmetric
            sData[idx11] = s3;

        } else {
            // Off-diagonal 4×4 system
            double[] A = lyapA4;
            double[] b = lyapB4;
            double[] x = lyapX4;

            // A is row-major 4×4.
            // Row 0
            A[0]  = d0 + e0;  A[1]  = e2;        A[2]  = d1;        A[3]  = 0.0;
            // Row 1
            A[4]  = e1;       A[5]  = d0 + e3;   A[6]  = 0.0;       A[7]  = d1;
            // Row 2
            A[8]  = d2;       A[9]  = 0.0;       A[10] = d3 + e0;   A[11] = e2;
            // Row 3
            A[12] = 0.0;      A[13] = d2;        A[14] = e1;        A[15] = d3 + e3;

            b[0] = v0;
            b[1] = v1;
            b[2] = v2;
            b[3] = v3;

            solve4x4InPlace(A, b, x);

            sData[idx00] = x[0];
            sData[idx01] = x[1];
            sData[idx10] = x[2];
            sData[idx11] = x[3];
        }
    }

    // ----------------------------------------------------------------------
    // Small linear solvers (3×3, 4×4)
    // ----------------------------------------------------------------------

    /**
     * Solve 3×3 system A x = b using Cramer's rule.
     *
     * A is flat row-major 3×3 of length 9.
     * b is length 3.
     * xOut is length 3 (overwritten).
     */
    private static void solve3x3InPlace(double[] A, double[] b, double[] xOut) {
        double a00 = A[0], a01 = A[1], a02 = A[2];
        double a10 = A[3], a11 = A[4], a12 = A[5];
        double a20 = A[6], a21 = A[7], a22 = A[8];

        double b0 = b[0], b1 = b[1], b2 = b[2];

        double det =
                a00 * (a11 * a22 - a12 * a21)
                        - a01 * (a10 * a22 - a12 * a20)
                        + a02 * (a10 * a21 - a11 * a20);

        if (Math.abs(det) < EPS_BLOCK) {
            throw new RuntimeException("Singular 3×3 Lyapunov system");
        }

        double det0 =
                b0  * (a11 * a22 - a12 * a21)
                        - a01 * (b1  * a22 - a12 * b2)
                        + a02 * (b1  * a21 - a11 * b2);

        double det1 =
                a00 * (b1  * a22 - a12 * b2)
                        - b0  * (a10 * a22 - a12 * a20)
                        + a02 * (a10 * b2  - b1  * a20);

        double det2 =
                a00 * (a11 * b2  - b1  * a21)
                        - a01 * (a10 * b2  - b1  * a20)
                        + b0  * (a10 * a21 - a11 * a20);

        double invDet = 1.0 / det;

        xOut[0] = det0 * invDet;
        xOut[1] = det1 * invDet;
        xOut[2] = det2 * invDet;
    }

    /**
     * Solve 4×4 system A x = b with Gaussian elimination + partial pivoting.
     *
     * A is flat row-major 4×4 of length 16.
     * b is length 4.
     * xOut is length 4 (overwritten).
     *
     * A and b are modified in-place during the solve.
     */
    private static void solve4x4InPlace(double[] A, double[] b, double[] xOut) {

        // Forward elimination
        for (int col = 0; col < 4; col++) {

            // Pivot selection
            int pivotRow = col;
            double maxVal = Math.abs(A[col * 4 + col]);
            for (int row = col + 1; row < 4; row++) {
                double val = Math.abs(A[row * 4 + col]);
                if (val > maxVal) {
                    maxVal = val;
                    pivotRow = row;
                }
            }

            if (maxVal < EPS_BLOCK) {
                throw new IllegalStateException("Singular 4×4 Lyapunov block");
            }

            // Row swap if needed
            if (pivotRow != col) {
                int baseCol   = col * 4;
                int basePivot = pivotRow * 4;
                for (int k = 0; k < 4; k++) {
                    double tmp = A[baseCol + k];
                    A[baseCol + k]   = A[basePivot + k];
                    A[basePivot + k] = tmp;
                }
                double tmpB = b[col];
                b[col]        = b[pivotRow];
                b[pivotRow]   = tmpB;
            }

            // Eliminate
            double pivot = A[col * 4 + col];
            for (int row = col + 1; row < 4; row++) {
                int baseRow = row * 4;
                double factor = A[baseRow + col] / pivot;
                if (factor == 0.0) {
                    continue;
                }
                for (int k = col; k < 4; k++) {
                    A[baseRow + k] -= factor * A[col * 4 + k];
                }
                b[row] -= factor * b[col];
            }
        }

        // Back substitution
        for (int i = 3; i >= 0; i--) {
            double sum = b[i];
            int baseRow = i * 4;
            for (int j = i + 1; j < 4; j++) {
                sum -= A[baseRow + j] * xOut[j];
            }
            double diag = A[baseRow + i];
            if (Math.abs(diag) < EPS_BLOCK) {
                throw new RuntimeException("Singular 4×4 Lyapunov block in back substitution");
            }
            xOut[i] = sum / diag;
        }
    }
}

