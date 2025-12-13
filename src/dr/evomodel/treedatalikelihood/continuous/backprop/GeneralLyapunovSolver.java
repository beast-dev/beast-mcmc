package dr.evomodel.treedatalikelihood.continuous.backprop;

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * Utility for solving Lyapunov equations of the form:
 *   S·X + X·S^T = RHS
 *
 * Uses Kronecker product formulation:
 *   (I ⊗ S + S ⊗ I) · vec(X) = vec(RHS)
 */
public class GeneralLyapunovSolver {
    private static final double EPS_BLOCK = 1e-12;

    /**
     * Solve S·X + X·S^T = RHS for X
     *
     * @param S_flat  flattened S matrix (row-major)
     * @param rhs     right-hand side matrix
     * @param dim     dimension
     * @return solution matrix X
     */
    public static DenseMatrix64F solve(double[] S_flat, DenseMatrix64F rhs, int dim) {

        // Build Kronecker matrix: (I ⊗ S + S ⊗ I)
        DenseMatrix64F A_L = new DenseMatrix64F(dim * dim, dim * dim);

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                int row_base = i * dim + j;

                for (int k = 0; k < dim; k++) {
                    for (int l = 0; l < dim; l++) {
                        int col = k * dim + l;

                        // (I ⊗ S)[row, col] = δ_{jl} · S[i,k]
                        double term1 = (j == l) ? S_flat[i * dim + k] : 0.0;

                        // (S ⊗ I)[row, col] = S[j,l] · δ_{ik}
                        double term2 = (i == k) ? S_flat[j * dim + l] : 0.0;

                        A_L.set(row_base, col, term1 + term2);
                    }
                }
            }
        }

        // Vectorize RHS
        DenseMatrix64F rhs_vec = new DenseMatrix64F(dim * dim, 1);
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                rhs_vec.set(i * dim + j, 0, rhs.get(i, j));
            }
        }

        // Solve linear system
        DenseMatrix64F solution = new DenseMatrix64F(dim * dim, 1);
        CommonOps.solve(A_L, rhs_vec, solution);

        // Reshape back to matrix
        DenseMatrix64F X = new DenseMatrix64F(dim, dim);
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                X.set(i, j, solution.get(i * dim + j, 0));
            }
        }

        return X;
    }

    /**
     * Compute dΣ_stat/dS[i,j] where S·Σ_stat + Σ_stat·S^T = Σ
     *
     * Solves: S·X + X·S^T = -(E_{ij}·Σ_stat + Σ_stat·E_{ij}^T)
     */
    public static DenseMatrix64F derivativeWrtS(double[] S_flat,
                                                DenseMatrix64F sigmaStat,
                                                int i_S, int j_S,
                                                int dim) {
        DenseMatrix64F rhs = new DenseMatrix64F(dim, dim);
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                double term1 = (i == i_S) ? sigmaStat.get(j_S, j) : 0.0;
                double term2 = (j == i_S) ? sigmaStat.get(i, j_S) : 0.0;
                rhs.set(i, j, -(term1 + term2));
            }
        }
        return solve(S_flat, rhs, dim);
    }

    /**
     * Compute dΣ_stat/dΣ[i,j] where S·Σ_stat + Σ_stat·S^T = Σ
     *
     * Solves: S·X + X·S^T = E_{ij}
     */
    public static DenseMatrix64F derivativeWrtSigma(double[] S_flat,
                                                    int i_Sig, int j_Sig,
                                                    int dim) {
        DenseMatrix64F rhs = new DenseMatrix64F(dim, dim);
        rhs.set(i_Sig, j_Sig, 1.0);
        return solve(S_flat, rhs, dim);
    }


    /**
     * Solve a single block equation: Dᵢ Σᵢⱼ + Σᵢⱼ Dⱼᵀ = Vᵢⱼ
     *
     * Cases:
     * - 1×1 × 1×1: σ = v / (dᵢ + dⱼ)
     * - 1×1 × 2×2: Solve 2×2 system
     * - 2×2 × 1×1: Solve 2×2 system
     * - 2×2 × 2×2: Solve 4×4 system analytically
     */
    public static void solveBlockLyapunov(DenseMatrix64F D,
                                          DenseMatrix64F V,
                                          DenseMatrix64F Sigma,
                                          int iStart, int iSize,
                                          int jStart, int jSize) {

        if (iSize == 1 && jSize == 1) {
            // 1×1 × 1×1 case: scalar equation
            // d_i * σ + σ * d_j = v  =>  σ = v / (d_i + d_j)
            double di = D.get(iStart, iStart);
            double dj = D.get(jStart, jStart);
            double v = V.get(iStart, jStart);

            double sigma = v / (di + dj);
            Sigma.set(iStart, jStart, sigma);

        } else if (iSize == 1 && jSize == 2) {
            // 1×1 × 2×2 case: 2 equations for 2 unknowns
            solve1x2BlockLyapunov(D, V, Sigma, iStart, jStart);

        } else if (iSize == 2 && jSize == 1) {
            // 2×1 × 1×1 case: 2 equations for 2 unknowns
            solve2x1BlockLyapunov(D, V, Sigma, iStart, jStart);

        } else if (iSize == 2 && jSize == 2) {
            // 2×2 × 2×2 case: 4 equations for 4 unknowns
            solve2x2BlockLyapunov(D, V, Sigma, iStart, jStart);
        }
    }

    /**
     * Solve: d * [σ₀, σ₁] + [σ₀, σ₁] * [d₀  d₁] = [v₀, v₁]
     *                                    [d₂  d₃]
     *
     * Expands to:
     *   (d + d₀) * σ₀ + d₂ * σ₁ = v₀
     *   d₁ * σ₀ + (d + d₃) * σ₁ = v₁
     */
    private static void solve1x2BlockLyapunov(DenseMatrix64F D,
                                              DenseMatrix64F V,
                                              DenseMatrix64F Sigma,
                                              int i, int j) {
        double d = D.get(i, i);
        double d0 = D.get(j, j);
        double d1 = D.get(j, j + 1);
        double d2 = D.get(j + 1, j);
        double d3 = D.get(j + 1, j + 1);

        double v0 = V.get(i, j);
        double v1 = V.get(i, j + 1);

        // Solve 2×2 system using Cramer's rule
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

        Sigma.set(i, j, s0);
        Sigma.set(i, j + 1, s1);
    }

    /**
     * Solve: [d₀  d₁] * [σ₀] + [σ₀] * d = [v₀]
     *        [d₂  d₃]   [σ₁]   [σ₁]       [v₁]
     *
     * Expands to:
     *   (d₀ + d) * σ₀ + d₁ * σ₁ = v₀
     *   d₂ * σ₀ + (d₃ + d) * σ₁ = v₁
     */
    private static void solve2x1BlockLyapunov(DenseMatrix64F D,
                                              DenseMatrix64F V,
                                              DenseMatrix64F Sigma,
                                              int i, int j) {
        double d0 = D.get(i, i);
        double d1 = D.get(i, i + 1);
        double d2 = D.get(i + 1, i);
        double d3 = D.get(i + 1, i + 1);
        double d = D.get(j, j);

        double v0 = V.get(i, j);
        double v1 = V.get(i + 1, j);

        // Solve 2×2 system using Cramer's rule
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

        Sigma.set(i, j, s0);
        Sigma.set(i + 1, j, s1);
    }

    /**
     * Solve: [d₀  d₁] * [σ₀  σ₁] + [σ₀  σ₁] * [e₀  e₁] = [v₀  v₁]
     *        [d₂  d₃]   [σ₂  σ₃]   [σ₂  σ₃]   [e₂  e₃]   [v₂  v₃]
     *
     * Special case: when i==j (diagonal block), Σ must be symmetric (σ₁ = σ₂),
     * reducing to a 3×3 system.
     */
    private static void solve2x2BlockLyapunov(DenseMatrix64F D,
                                              DenseMatrix64F V,
                                              DenseMatrix64F Sigma,
                                              int i, int j) {
        // Extract Di block
        double d0 = D.get(i, i);
        double d1 = D.get(i, i + 1);
        double d2 = D.get(i + 1, i);
        double d3 = D.get(i + 1, i + 1);

        // Extract Dj block (transposed in equation)
        double e0 = D.get(j, j);
        double e1 = D.get(j, j + 1);
        double e2 = D.get(j + 1, j);
        double e3 = D.get(j + 1, j + 1);

        // Extract V block
        double v0 = V.get(i, j);
        double v1 = V.get(i, j + 1);
        double v2 = V.get(i + 1, j);
        double v3 = V.get(i + 1, j + 1);

        if (i == j) {
            // Diagonal block: Σ must be symmetric
            // D * Σ + Σ * Dᵀ = V with Σ = Σᵀ
            //
            // [d₀  d₁] [σ₀  σ₁]   [σ₀  σ₁] [d₀  d₂]   [v₀  v₁]
            // [d₂  d₃] [σ₁  σ₃] + [σ₁  σ₃] [d₁  d₃] = [v₂  v₃]
            //
            // (0,0): d₀σ₀ + d₁σ₁ + σ₀d₀ + σ₁d₁ = v₀  =>  2d₀σ₀ + 2d₁σ₁ = v₀
            // (0,1): d₀σ₁ + d₁σ₃ + σ₀d₂ + σ₁d₃ = v₁  =>  d₂σ₀ + (d₀+d₃)σ₁ + d₁σ₃ = v₁
            // (1,1): d₂σ₁ + d₃σ₃ + σ₁d₂ + σ₃d₃ = v₃  =>  2d₂σ₁ + 2d₃σ₃ = v₃

            double[][] A_mat = new double[3][3];
            A_mat[0][0] = 2*d0;        A_mat[0][1] = 2*d1;      A_mat[0][2] = 0;
            A_mat[1][0] = d2;          A_mat[1][1] = d0+d3;     A_mat[1][2] = d1;
            A_mat[2][0] = 0;           A_mat[2][1] = 2*d2;      A_mat[2][2] = 2*d3;

            double[] b_vec = {v0, v1, v3};
            double[] x = solve3x3(A_mat, b_vec);

            double s0 = x[0];
            double s1 = x[1];
            double s3 = x[2];

            Sigma.set(i,     j,     s0);
            Sigma.set(i,     j + 1, s1);
            Sigma.set(i + 1, j,     s1);  // Symmetric
            Sigma.set(i + 1, j + 1, s3);

        } else {
            // Off-diagonal block: no symmetry constraint
            // Solve full 4×4 system
            //
            // Row 0 (eq for V[0,0]): (d0+e0)*σ₀ + e2*σ₁ + d1*σ₂ + 0*σ₃ = v₀
            // Row 1 (eq for V[0,1]): e1*σ₀ + (d0+e3)*σ₁ + 0*σ₂ + d1*σ₃ = v₁
            // Row 2 (eq for V[1,0]): d2*σ₀ + 0*σ₁ + (d3+e0)*σ₂ + e2*σ₃ = v₂
            // Row 3 (eq for V[1,1]): 0*σ₀ + d2*σ₁ + e1*σ₂ + (d3+e3)*σ₃ = v₃

            double[][] A = new double[4][4];
            A[0][0] = d0 + e0; A[0][1] = e2;      A[0][2] = d1;      A[0][3] = 0;
            A[1][0] = e1;      A[1][1] = d0 + e3; A[1][2] = 0;       A[1][3] = d1;
            A[2][0] = d2;      A[2][1] = 0;       A[2][2] = d3 + e0; A[2][3] = e2;
            A[3][0] = 0;       A[3][1] = d2;      A[3][2] = e1;      A[3][3] = d3 + e3;

            double[] b = {v0, v1, v2, v3};
            double[] x = solve4x4(A, b);

            Sigma.set(i,     j,     x[0]);
            Sigma.set(i,     j + 1, x[1]);
            Sigma.set(i + 1, j,     x[2]);
            Sigma.set(i + 1, j + 1, x[3]);
        }
    }

    /**
     * Solve 3×3 system using Cramer's rule
     */
    private static double[] solve3x3(double[][] A, double[] b) {
        // Compute determinant
        double det = A[0][0] * (A[1][1]*A[2][2] - A[1][2]*A[2][1])
                - A[0][1] * (A[1][0]*A[2][2] - A[1][2]*A[2][0])
                + A[0][2] * (A[1][0]*A[2][1] - A[1][1]*A[2][0]);

        if (Math.abs(det) < EPS_BLOCK) {
            throw new RuntimeException("Singular 3×3 Lyapunov system");
        }

        // Cramer's rule
        double[] x = new double[3];

        // x[0]
        double det0 = b[0] * (A[1][1]*A[2][2] - A[1][2]*A[2][1])
                - A[0][1] * (b[1]*A[2][2] - A[1][2]*b[2])
                + A[0][2] * (b[1]*A[2][1] - A[1][1]*b[2]);
        x[0] = det0 / det;

        // x[1]
        double det1 = A[0][0] * (b[1]*A[2][2] - A[1][2]*b[2])
                - b[0] * (A[1][0]*A[2][2] - A[1][2]*A[2][0])
                + A[0][2] * (A[1][0]*b[2] - b[1]*A[2][0]);
        x[1] = det1 / det;

        // x[2]
        double det2 = A[0][0] * (A[1][1]*b[2] - b[1]*A[2][1])
                - A[0][1] * (A[1][0]*b[2] - b[1]*A[2][0])
                + b[0] * (A[1][0]*A[2][1] - A[1][1]*A[2][0]);
        x[2] = det2 / det;

        return x;
    }

    /**
     * Solve 4×4 system using Gaussian elimination with partial pivoting.
     * Optimized for small fixed size.
     */
    private static double[] solve4x4(double[][] A, double[] b) {
        // Copy to avoid modifying input
        double[][] M = new double[4][4];
        double[] rhs = new double[4];
        for (int i = 0; i < 4; i++) {
            System.arraycopy(A[i], 0, M[i], 0, 4);
            rhs[i] = b[i];
        }

        // Gaussian elimination with partial pivoting
        for (int col = 0; col < 4; col++) {
            // Find pivot
            int pivotRow = col;
            double maxVal = Math.abs(M[col][col]);
            for (int row = col + 1; row < 4; row++) {
                double val = Math.abs(M[row][col]);
                if (val > maxVal) {
                    maxVal = val;
                    pivotRow = row;
                }
            }

            if (maxVal < EPS_BLOCK) {
                throw new RuntimeException("Singular 4×4 Lyapunov block");
            }

            // Swap rows if needed
            if (pivotRow != col) {
                double[] tempRow = M[col];
                M[col] = M[pivotRow];
                M[pivotRow] = tempRow;

                double tempRhs = rhs[col];
                rhs[col] = rhs[pivotRow];
                rhs[pivotRow] = tempRhs;
            }

            // Eliminate
            for (int row = col + 1; row < 4; row++) {
                double factor = M[row][col] / M[col][col];
                for (int k = col; k < 4; k++) {
                    M[row][k] -= factor * M[col][k];
                }
                rhs[row] -= factor * rhs[col];
            }
        }

        // Back substitution
        double[] x = new double[4];
        for (int i = 3; i >= 0; i--) {
            double sum = rhs[i];
            for (int j = i + 1; j < 4; j++) {
                sum -= M[i][j] * x[j];
            }
            x[i] = sum / M[i][i];
        }

        return x;
    }
}