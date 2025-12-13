package dr.evomodel.treedatalikelihood.continuous.backprop;

import org.ejml.data.DenseMatrix64F;

public class BlockDiagonalExpSolver {
    static double EPS_BLOCK = 1e-12;
    /**
     * Compute Q_eig(t) = exp(-D t) for a block-tridiagonal D whose blocks
     * are 1×1 and general 2×2.
//     *
     * Works directly with blockDParams array, avoiding D matrix construction.
     *
     * Contiguous nonzero off-diagonal pairs (i, i+1) define 2×2 blocks:
     *   B_i = [ d_i      u_i     ]
     *         [ l_i   d_{i+1}    ].
     * Otherwise we treat d_i as a 1×1 block.
     */
    public static void expBlockMatrix(double[] blockDParams, int dim, double t,
                                      DenseMatrix64F Q) {

        final int upperOffset = dim;
        final int lowerOffset = dim + (dim - 1);

        double[] qData = Q.data;
        java.util.Arrays.fill(qData, 0.0);

        int i = 0;
        while (i < dim) {
            double u = (i < dim - 1) ? blockDParams[upperOffset + i] : 0.0;
            double l = (i < dim - 1) ? blockDParams[lowerOffset + i] : 0.0;

            if (i == dim - 1 || (Math.abs(u) < EPS_BLOCK && Math.abs(l) < EPS_BLOCK)) {
                // 1×1 block
                double a = blockDParams[i];
                int idx = i * dim + i;
                qData[idx] = Math.exp(-a * t);
                i += 1;
            } else {
                // General 2×2 block:
                //   [ a  b ]
                //   [ c  d ]
                double a = blockDParams[i];
                double d = blockDParams[i + 1];
                double b = u;
                double c = l;

                fillExp2x2GeneralBlock(a, b, c, d, t, qData, dim, i);
                i += 2;
            }
        }
    }

    /**
     * Fill the 2×2 block of Q corresponding to exp(-t B), where
     *
     *   B = [ a  b ]
     *       [ c  d ].
     *
     * We use the standard 2×2 formula based on
     *   τ  = (a + d)/2,
     *   Δ² = ((a - d)/2)² + b c,
     *   K  = B - τ I (tr(K) = 0, K² = Δ² I),
     *
     * and
     *   exp(-t B) = e^{-τ t} [ α I + β K ],
     *
     * with
     *   (Δ² > 0): α = cosh(δ t), β = -sinh(δ t)/δ
     *   (Δ² < 0): α = cos(δ t),  β = -sin(δ t)/δ
     *   (Δ² = 0): α = 1,         β = -t.
     *
     * Writes directly into the flat qData array in row-major order.
     */
    private static void fillExp2x2GeneralBlock(double a,
                                               double b,
                                               double c,
                                               double d,
                                               double t,
                                               double[] qData,
                                               int dim,
                                               int offsetRowCol) {

        double tau    = 0.5 * (a + d);
        double diff   = 0.5 * (a - d);
        double delta2 = diff * diff + b * c;

        double alpha;
        double beta;

        if (Math.abs(delta2) < EPS_BLOCK) {
            // Δ² ≈ 0  ⇒  K² ≈ 0 and exp(-t B) ≈ e^{-τ t}(I - t K)
            alpha = 1.0;
            beta  = -t;
        } else if (delta2 > 0.0) {
            double delta = Math.sqrt(delta2);
            double dtv   = delta * t;
            alpha = Math.cosh(dtv);
            beta  = -Math.sinh(dtv) / delta;
        } else {
            double delta = Math.sqrt(-delta2);
            double dtv   = delta * t;
            alpha = Math.cos(dtv);
            beta  = -Math.sin(dtv) / delta;
        }

        double factor = Math.exp(-tau * t);

        // K = B - τ I
        double k11 = a - tau;
        double k22 = d - tau;

        double m11 = alpha + beta * k11;
        double m12 = beta * b;
        double m21 = beta * c;
        double m22 = alpha + beta * k22;

        int i = offsetRowCol;
        int j = offsetRowCol + 1;

        int idx00 = i       * dim + i;
        int idx01 = i       * dim + j;
        int idx10 = (i + 1) * dim + i;
        int idx11 = (i + 1) * dim + j;

        qData[idx00] = factor * m11;
        qData[idx01] = factor * m12;
        qData[idx10] = factor * m21;
        qData[idx11] = factor * m22;
    }
}
