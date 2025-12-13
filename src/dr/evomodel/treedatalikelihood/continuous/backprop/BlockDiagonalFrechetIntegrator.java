package dr.evomodel.treedatalikelihood.continuous.backprop;

/**
 * Low-level helper for 2×2-based Fréchet integrals on block-diagonalizable OU models.
 *
 * All methods are allocation-free: results are written into caller-provided buffers.
 *
 * Expected buffer sizes:
 *  - out2: length >= 2  (for 1×2 or 2×1 results)
 *  - out4 / expBuf: length >= 4 (for 2×2 matrices in row-major order)
 *
 * @author Filippo Monti
 */
public final class BlockDiagonalFrechetIntegrator {

    // 5-point Gauss–Legendre on [0,1] (nodes and weights)
    private static final double[] GL5_NODES = {
            0.046910077030668,
            0.230765344947158,
            0.500000000000000,
            0.769234655052842,
            0.953089922969332
    };

    private static final double[] GL5_WEIGHTS = {
            0.118463442528095,
            0.239314335249683,
            0.284444444444444,
            0.239314335249683,
            0.118463442528095
    };

    private static final double SING_EPS = 1e-10;

    private BlockDiagonalFrechetIntegrator() {
        // no instances
    }

    /**
     * Compute exponential of 2×2 matrix [[a, b], [c, d]].
     * Result is written into {@code out} in row-major order:
     *   out[0] = exp00, out[1] = exp01,
     *   out[2] = exp10, out[3] = exp11.
     *
     * Uses the analytical formula for the 2×2 matrix exponential.
     */
    public static void exp2x2(double a, double b, double c, double d,
                              double[] out) {
        double tau = 0.5 * (a + d);
        double diff = 0.5 * (a - d);
        double delta2 = diff * diff + b * c;

        double alpha, beta;

        if (Math.abs(delta2) < 1e-12) {
            // Δ² ≈ 0
            alpha = 1.0;
            beta = 1.0;
        } else if (delta2 > 0.0) {
            double delta = Math.sqrt(delta2);
            alpha = Math.cosh(delta);
            beta = Math.sinh(delta) / delta;
        } else {
            double delta = Math.sqrt(-delta2);
            alpha = Math.cos(delta);
            beta = Math.sin(delta) / delta;
        }

        double factor = Math.exp(tau);

        // K = [[a-τ, b], [c, d-τ]]
        double k11 = a - tau;
        double k22 = d - tau;

        // exp = factor * (alpha*I + beta*K)
        out[0] = factor * (alpha + beta * k11); // exp00
        out[1] = factor * (beta * b);           // exp01
        out[2] = factor * (beta * c);           // exp10
        out[3] = factor * (alpha + beta * k22); // exp11
    }

    /**
     * Integrate ∫₀¹ exp(sλ_i) [e1, e2] exp((1-s)A_j) ds
     * where A_j = [[a, b], [c, d]] is 2×2.
     *
     * Exact formula via 2×2 linear system:
     *
     *   M = -A_j^T + λ_i I₂ = [[λ_i - a,   -c   ],
     *                          [  -b,    λ_i - d]]
     *
     * Let expM = exp(M), E^T = [e1; e2].
     * Then RHS = (expM - I₂) E^T, and y solves M y = RHS.
     * The integral is then (exp(A_j^T) y)^T.
     *
     * If |det(M)| < SING_EPS, falls back to 5-point Gauss–Legendre quadrature.
     *
     * @param lambda_i scalar eigenvalue (already scaled by -t)
     * @param a,b,c,d  entries of A_j (already scaled by -t)
     * @param e1,e2    entries of 1×2 block E_{ij}
     * @param out      length ≥ 2, will hold the resulting 1×2 row
     * @param expBuf   length ≥ 4, scratch for 2×2 exponentials
     */
    public static void integrate1x1_2x2(double lambda_i,
                                        double a, double b, double c, double d,
                                        double e1, double e2,
                                        double[] out,
                                        double[] expBuf) {
        // Reset output
        out[0] = 0.0;
        out[1] = 0.0;

        // M = -A_j^T + λ_i I₂ = [[λ_i - a, -c],
        //                        [  -b,   λ_i - d]]
        double m00 = lambda_i - a;
        double m01 = -c;
        double m10 = -b;
        double m11 = lambda_i - d;

        double detM = m00 * m11 - m01 * m10;

        if (Math.abs(detM) < SING_EPS) {
            // Fallback: Gauss–Legendre quadrature (no allocations)
            for (int q = 0; q < GL5_NODES.length; q++) {
                double s = GL5_NODES[q];
                double w = GL5_WEIGHTS[q];

                double exp_s_lambda_i = Math.exp(s * lambda_i);

                // exp((1-s) A_j)
                double scale = (1.0 - s);
                exp2x2(scale * a, scale * b, scale * c, scale * d, expBuf);

                // [e1, e2] * exp((1-s)A_j)
                double prod1 = e1 * expBuf[0] + e2 * expBuf[2];
                double prod2 = e1 * expBuf[1] + e2 * expBuf[3];

                out[0] += w * exp_s_lambda_i * prod1;
                out[1] += w * exp_s_lambda_i * prod2;
            }
            return;
        }

        // exp(M)
        exp2x2(m00, m01, m10, m11, expBuf);

        // (exp(M) - I) * E^T
        double rhs0 = (expBuf[0] - 1.0) * e1 + expBuf[1] * e2;
        double rhs1 = expBuf[2] * e1 + (expBuf[3] - 1.0) * e2;

        // Solve M y = RHS (2×2 analytic inverse)
        double invDet = 1.0 / detM;
        double y0 = ( rhs0 * m11 - rhs1 * m01) * invDet;
        double y1 = (-rhs0 * m10 + rhs1 * m00) * invDet;

        // exp(A_j^T)
        exp2x2(a, c, b, d, expBuf);

        // integral_col = exp(A_j^T) * y
        out[0] = expBuf[0] * y0 + expBuf[1] * y1;
        out[1] = expBuf[2] * y0 + expBuf[3] * y1;
    }

    /**
     * Integrate ∫₀¹ exp(sA_i) [e1; e2] exp((1-s)λ_j) ds
     * where A_i = [[a, b], [c, d]] and λ_j is scalar.
     *
     * Exact formula:
     *
     *   M = A_i - λ_j I₂ = [[a - λ_j, b],
     *                       [c,       d - λ_j]]
     *
     *   exp(sA_i) E exp((1-s)λ_j)
     *     = exp(λ_j) exp(sM) E
     *
     * ⇒ integral = exp(λ_j) (exp(M) - I) M^{-1} E.
     *
     * If |det(M)| < SING_EPS, falls back to 5-point Gauss–Legendre.
     *
     * @param a,b,c,d   entries of A_i (already scaled by -t)
     * @param lambda_j  scalar eigenvalue (already scaled by -t)
     * @param e1,e2     entries of 2×1 block E_{ij}
     * @param out       length ≥ 2, will hold the resulting 2×1 column
     * @param expBuf    length ≥ 4, scratch for 2×2 exponentials
     */
    public static void integrate2x2_1x1(double a, double b, double c, double d,
                                        double lambda_j,
                                        double e1, double e2,
                                        double[] out,
                                        double[] expBuf) {
        // Reset output
        out[0] = 0.0;
        out[1] = 0.0;

        // M = A_i - λ_j I₂
        double m00 = a - lambda_j;
        double m01 = b;
        double m10 = c;
        double m11 = d - lambda_j;

        double detM = m00 * m11 - m01 * m10;
        if (Math.abs(detM) < SING_EPS) {
            // Fallback: Gauss–Legendre
            for (int q = 0; q < GL5_NODES.length; q++) {
                double s = GL5_NODES[q];
                double w = GL5_WEIGHTS[q];

                // exp(s A_i)
                double scale = s;
                exp2x2(scale * a, scale * b, scale * c, scale * d, expBuf);
                double exp_1ms_lambda_j = Math.exp((1.0 - s) * lambda_j);

                double prod1 = expBuf[0] * e1 + expBuf[1] * e2;
                double prod2 = expBuf[2] * e1 + expBuf[3] * e2;

                out[0] += w * prod1 * exp_1ms_lambda_j;
                out[1] += w * prod2 * exp_1ms_lambda_j;
            }
            return;
        }

        // exp(M)
        exp2x2(m00, m01, m10, m11, expBuf);

        // (exp(M) - I) * E
        double rhs0 = (expBuf[0] - 1.0) * e1 + expBuf[1] * e2;
        double rhs1 = expBuf[2] * e1 + (expBuf[3] - 1.0) * e2;

        // Solve M y = RHS
        double invDet = 1.0 / detM;
        double y0 = ( rhs0 * m11 - rhs1 * m01) * invDet;
        double y1 = (-rhs0 * m10 + rhs1 * m00) * invDet;

        double factor = Math.exp(lambda_j);
        out[0] = factor * y0;
        out[1] = factor * y1;
    }

    /**
     * Integrate ∫₀¹ exp(sA_i) E exp((1-s)A_j) ds
     * where both A_i and A_j are 2×2 matrices:
     *
     *   A_i = [[ai, bi],
     *          [ci, di]]
     *
     *   A_j = [[aj, bj],
     *          [cj, dj]]
     *
     * E is 2×2 with entries (e11, e12; e21, e22).
     *
     * This implementation uses 5-point Gauss–Legendre quadrature only.
     * For 2×2 blocks this is typically competitive or faster than a 4×4
     * Kronecker-sum exact path, and is significantly simpler.
     *
     * @param ai,bi,ci,di entries of A_i (already scaled by -t)
     * @param aj,bj,cj,dj entries of A_j (already scaled by -t)
     * @param e11,e12,e21,e22 entries of E_{ij}
     * @param out   length ≥ 4, will hold 2×2 result in row-major order
     * @param expAi length ≥ 4, scratch for exp(s A_i)
     * @param expAj length ≥ 4, scratch for exp((1-s) A_j)
     */
    public static void integrate2x2_2x2(double ai, double bi, double ci, double di,
                                        double aj, double bj, double cj, double dj,
                                        double e11, double e12, double e21, double e22,
                                        double[] out,
                                        double[] expAi,
                                        double[] expAj) {
        // Reset output
        out[0] = out[1] = out[2] = out[3] = 0.0;

        for (int q = 0; q < GL5_NODES.length; q++) {
            double s = GL5_NODES[q];
            double w = GL5_WEIGHTS[q];

            // exp(s A_i)
            double scale_i = s;
            exp2x2(scale_i * ai, scale_i * bi, scale_i * ci, scale_i * di, expAi);

            // exp((1-s) A_j)
            double scale_j = 1.0 - s;
            exp2x2(scale_j * aj, scale_j * bj, scale_j * cj, scale_j * dj, expAj);

            // T = exp(s A_i) * E
            double t00 = expAi[0] * e11 + expAi[1] * e21;
            double t01 = expAi[0] * e12 + expAi[1] * e22;
            double t10 = expAi[2] * e11 + expAi[3] * e21;
            double t11 = expAi[2] * e12 + expAi[3] * e22;

            // R = T * exp((1-s) A_j)
            double r00 = t00 * expAj[0] + t01 * expAj[2];
            double r01 = t00 * expAj[1] + t01 * expAj[3];
            double r10 = t10 * expAj[0] + t11 * expAj[2];
            double r11 = t10 * expAj[1] + t11 * expAj[3];

            out[0] += w * r00;
            out[1] += w * r01;
            out[2] += w * r10;
            out[3] += w * r11;
        }
    }
}
