package dr.evomodel.treedatalikelihood.continuous.backprop;

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * Helper for the adjoint Lyapunov contribution in OU backprop when
 * S = R D R^{-1} with D block diagonal.
 *
 * It works in the D-basis, using BlockDiagonalLyapunovSolver on the
 * compressed representation of D^T, then maps back to the S-basis
 * and accumulates ∇_S L in the caller's dLdS.
 *
 * NOT thread-safe: owns internal k×k workspaces.
 */

/**
 * @author Filippo Monti
 */
final class BlockDiagonalLyapunovAdjointHelper {

    private final int k;
    private final BlockDiagonalLyapunovSolver lyapSolver;

    // R and R^{-1}
    private final DenseMatrix64F R;
    private final DenseMatrix64F Rinv;

    // Compressed representation of D^T:
    // [d_0..d_{k-1}, u_0..u_{k-2}, l_0..l_{k-2}]
    private final double[] blockDAdjointParams;

    // Workspaces (all k×k)
    private final DenseMatrix64F tmpK1;
    private final DenseMatrix64F tmpK2;
    private final DenseMatrix64F tmpK3;
    private final DenseMatrix64F tmpK4;
    private final DenseMatrix64F tmpK5;

    BlockDiagonalLyapunovAdjointHelper(int k,
                                       BlockDiagonalLyapunovSolver lyapSolver,
                                       DenseMatrix64F R,
                                       DenseMatrix64F Rinv) {
        this.k = k;
        this.lyapSolver = lyapSolver;
        this.R = R;
        this.Rinv = Rinv;

        this.blockDAdjointParams = new double[3 * k - 2];

        this.tmpK1 = new DenseMatrix64F(k, k);
        this.tmpK2 = new DenseMatrix64F(k, k);
        this.tmpK3 = new DenseMatrix64F(k, k);
        this.tmpK4 = new DenseMatrix64F(k, k);
        this.tmpK5 = new DenseMatrix64F(k, k);
    }

    /**
     * Accumulate the Lyapunov adjoint contribution into dLdSAcc:
     *
     *   D^T · Ỹ + Ỹ · D = -G̃_X   in D-basis,
     *
     * where G̃_X = R^{-1} (dL/dΣ_stat) R^{-T},
     * Y = R^T Ỹ R^{-1},
     * and finally
     *
     *   ∇_S L += Y Σ_stat + Y^T Σ_stat^T. (recall Σ_stat symmetric)
     *
     * @param Sigma_stat      stationary covariance Σ_stat (S-basis)
     * @param dLdSigmaStat    gradient dL/dΣ_stat (S-basis)
     * @param blockDParams    compressed D params [diag | upper | lower]
     * @param dLdSAcc         accumulator for ∇_S L (S-basis)
     */
    void accumulateLyapunovContribution(DenseMatrix64F Sigma_stat,
                                        DenseMatrix64F dLdSigmaStat,
                                        double[] blockDParams,
                                        DenseMatrix64F dLdSAcc) {

        final int diagLen = k;
        final int offLen  = k - 1;

        // Build compressed D^T from D:
        //
        // D[i,i]   = d_i
        // D[i,i+1] = u_i
        // D[i+1,i] = l_i
        //
        // D^T[i,i]   = d_i
        // D^T[i,i+1] = l_i
        // D^T[i+1,i] = u_i
        //
        // => swap upper and lower bands
        System.arraycopy(blockDParams, 0,
                blockDAdjointParams, 0,
                diagLen);

        System.arraycopy(blockDParams, diagLen + offLen,
                blockDAdjointParams, diagLen,
                offLen); // upper(D^T) = lower(D)

        System.arraycopy(blockDParams, diagLen,
                blockDAdjointParams, diagLen + offLen,
                offLen); // lower(D^T) = upper(D)

        // G̃_X = R^{-1} · dLdSigmaStat · (Rinv)^T
        DenseMatrix64F temp          = tmpK1;
        DenseMatrix64F G_tilde_X     = tmpK2;
        DenseMatrix64F neg_G_tilde_X = tmpK3;
        DenseMatrix64F Y_tilde       = tmpK4;

        CommonOps.mult(Rinv, dLdSigmaStat, temp);
        CommonOps.multTransB(temp, Rinv, G_tilde_X);

        // V = -G̃_X
        CommonOps.scale(-1.0, G_tilde_X, neg_G_tilde_X);

        // Solve: D^T Ỹ + Ỹ D = V
        lyapSolver.solve(blockDAdjointParams, neg_G_tilde_X, Y_tilde);

        // Map back to S-basis: Y = R^T Ỹ R^{-1}
        DenseMatrix64F Y = tmpK5;
        CommonOps.multTransA(R, Y_tilde, temp); // temp = R^T Ỹ
        CommonOps.mult(temp, Rinv, Y);          // Y = (R^T Ỹ) R^{-1}

        DenseMatrix64F YT_Sigma = tmpK3;  // Y^T Σ

        CommonOps.multAdd(1.0, Y, Sigma_stat, dLdSAcc);   // += Y Σ
        CommonOps.multTransA(Y, Sigma_stat, YT_Sigma);   // Y^T Σ
        CommonOps.addEquals(dLdSAcc, YT_Sigma);           // += Y^T Σ



//        CommonOps.mult(Y,         Sigma_stat, term1); // term1 = Y Σ
//        CommonOps.mult(Sigma_stat, Y,         tmp);   // tmp   = Σ Y
//        CommonOps.transpose(tmp, term2);              // term2 = (Σ Y)^T
//
//        CommonOps.add(term1, term2, contrib_4);
//
//        // Accumulate into dLdS
//        CommonOps.addEquals(dLdSAcc, contrib_4);
    }
}
