package dr.evomodel.treedatalikelihood.continuous.backprop;

import dr.evomodel.continuous.MultivariateElasticModel;
import dr.inference.model.BlockDiagonalCosSinMatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * OU backprop strategy that computes primitive gradients:
 *   - dL/dS
 *   - dL/dSigma_stat
 *   - dL/dMu
 *
 * It does NOT map dL/dS into gradients w.r.t. (R, a, r, theta).
 * That mapping must is done by a separate mapper (e.g., BlockDiagCosSinPrimitiveGradientMapper).
 *
 * Assumes S = R D R^{-1}, with D block-diagonal (1x1 and 2x2 blocks) and general invertible R.
 *
 * Not thread-safe (reuses workspaces).
 */
/*
* @author Filippo Monti
 */
public final class BlockDiagonalizableBackpropStrategy implements PrimitiveParameterBackpropStrategy {

    private final BlockDiagonalCosSinMatrixParameter blockParam;

    // Current R and R^{-1} (refreshed each call)
    private final DenseMatrix64F R;
    private final DenseMatrix64F Rinv;

    private final int k;
    private final boolean debug;

    // Lyapunov adjoint helper (uses R, Rinv references)
    private final BlockDiagonalLyapunovAdjointHelper lyapAdjHelper;

    // Reusable k×k workspaces (not thread safe)
    private final DenseMatrix64F tmpK1;
    private final DenseMatrix64F tmpK2;
    private final DenseMatrix64F tmpK3;
    private final DenseMatrix64F tmpK4;

    // Reusable accumulators (avoid per-call allocations)
    private final DenseMatrix64F dLdS;
    private final DenseMatrix64F dLdMu;
    private final DenseMatrix64F dLdSigmaStat;

    // Reusable compressed D params: [diag (k), upper (k-1), lower (k-1)]
    private final double[] blockDParams;

    // Helper that handles Fréchet adjoint integrals in the D-basis (holds references to R, Rinv)
    private final BlockDiagonalFrechetHelper frechetHelper;

    public BlockDiagonalizableBackpropStrategy(MultivariateElasticModel elasticModel) {

        if (elasticModel == null) {
            throw new IllegalArgumentException("elasticModel cannot be null");
        }

        MatrixParameterInterface mpi = elasticModel.getMatrixParameterInterface();
        if (!(mpi instanceof BlockDiagonalCosSinMatrixParameter)) {
            throw new IllegalArgumentException(
                    "BlockDiagonalizableBackpropStrategy requires a BlockDiagonalCosSinMatrixParameter"
            );
        }
        this.blockParam = (BlockDiagonalCosSinMatrixParameter) mpi;

        this.k = blockParam.getRowDimension();
        this.debug = false;

        // R and R^{-1}
        this.R = new DenseMatrix64F(k, k);
        this.Rinv = new DenseMatrix64F(k, k);
        blockParam.fillRAndRinv(R.getData(), Rinv.getData());

        // Workspaces / accumulators
        this.tmpK1 = new DenseMatrix64F(k, k);
        this.tmpK2 = new DenseMatrix64F(k, k);
        this.tmpK3 = new DenseMatrix64F(k, k);
        this.tmpK4 = new DenseMatrix64F(k, k);

        this.dLdS = new DenseMatrix64F(k, k);
        this.dLdMu = new DenseMatrix64F(k, 1);
        this.dLdSigmaStat = new DenseMatrix64F(k, k);

        this.blockDParams = new double[3 * k - 2];

        // Lyapunov solver/helper
        BlockDiagonalLyapunovSolver solver = new BlockDiagonalLyapunovSolver(k);
        this.lyapAdjHelper = new BlockDiagonalLyapunovAdjointHelper(k, solver, R, Rinv);

        // Fréchet helper (all math lives there; it uses R and Rinv by reference)
        this.frechetHelper = new BlockDiagonalFrechetHelper(blockParam, R, Rinv);

        if (debug) {
            System.out.println("BlockDiagonalizableBackpropStrategy initialized:");
            System.out.println("  k = " + k);
        }
    }

    @Override
    public PrimitiveGradientSet backprop(
            ContinuousTraitBackpropGradient.OUBranchCache cache,
            MessageBackprop.Result leafGrads) {

        if (debug) {
            System.out.println("\n=== BlockDiagonalizableBackpropStrategy.backprop ===");
        }

        // Pull cached OU quantities
        final DenseMatrix64F Sigma_stat = cache.sigmaStat;
        final DenseMatrix64F A          = cache.A;
        final DenseMatrix64F mu         = cache.mu;
        final double         t          = cache.t;

        // Pull leaf gradients
        final DenseMatrix64F dLdA = leafGrads.dLdA;
        final DenseMatrix64F dLdV = leafGrads.dLdV;
        final DenseMatrix64F dLdb = leafGrads.dLdb;

        // Refresh R and R^{-1} (R changes in MCMC)
        blockParam.fillRAndRinv(R.getData(), Rinv.getData());

        // Refresh compressed D params without allocating
        // If you don't trust fillBlockDiagonalElements yet, replace with:
        // System.arraycopy(blockParam.getBlockDiagonalElements(), 0, blockDParams, 0, blockDParams.length);
        blockParam.fillBlockDiagonalElements(blockDParams);

        // Zero accumulators
        dLdS.zero();
        dLdMu.zero();
        dLdSigmaStat.zero();

        // 1) dL/dMu and dL/dSigmaStat from b and V bookkeeping
        computePreliminaryGradients(A, dLdb, dLdV, dLdMu, dLdSigmaStat);

        // 2) Contributions through exp(-St)
        computeExponentialContributions(Sigma_stat, A, mu, t, dLdA, dLdV, dLdb, dLdS, blockDParams);

        // 3) Contribution through Lyapunov (Sigma_stat depends on S)
        computeLyapunovContribution(Sigma_stat, dLdSigmaStat, dLdS, blockDParams);

        if (debug) {
            System.out.println("=== BlockDiagonalizableBackpropStrategy complete ===");
        }

        // IMPORTANT:
        // - dLdS is w.r.t. S (not w.r.t. parameterization).
        // - Mapping to (R, a, r, theta) must be done elsewhere.
        return new PrimitiveGradientSet(dLdS, null, dLdSigmaStat, dLdMu);
    }

    // ---------- 1. Preliminary gradients ----------

    private void computePreliminaryGradients(DenseMatrix64F A,
                                             DenseMatrix64F dLdb,
                                             DenseMatrix64F dLdV,
                                             DenseMatrix64F dLdMuOut,
                                             DenseMatrix64F dLdSigmaStatOut) {
        // I_minus_A = I - A
        DenseMatrix64F I_minus_A = tmpK1;
        CommonOps.setIdentity(I_minus_A);
        CommonOps.subtractEquals(I_minus_A, A);

        // dLdMu = (I-A)^T dLdb
        DenseMatrix64F I_minus_A_T = tmpK2;
        CommonOps.transpose(I_minus_A, I_minus_A_T);
        CommonOps.mult(I_minus_A_T, dLdb, dLdMuOut);

        // dLdSigmaStat = dLdV - A^T dLdV A
        DenseMatrix64F A_T = tmpK3;
        CommonOps.transpose(A, A_T);

        DenseMatrix64F temp = tmpK4;
        CommonOps.mult(A_T, dLdV, temp);
        CommonOps.mult(temp, A, dLdSigmaStatOut);
        CommonOps.scale(-1.0, dLdSigmaStatOut);
        CommonOps.addEquals(dLdSigmaStatOut, dLdV);
    }

    // ---------- 2. Contributions 1–3: Fréchet via helper ----------

    private void computeExponentialContributions(DenseMatrix64F Sigma_stat,
                                                 DenseMatrix64F A,
                                                 DenseMatrix64F mu,
                                                 double t,
                                                 DenseMatrix64F dLdA,
                                                 DenseMatrix64F dLdV,
                                                 DenseMatrix64F dLdb,
                                                 DenseMatrix64F dLdSAcc,
                                                 double[] blockDParams) {

        // 1) A = exp(-S t)
        DenseMatrix64F Z1 = frechetHelper.frechetAdjointExp(blockDParams, dLdA, t);
        CommonOps.scale(-t, Z1);
        CommonOps.addEquals(dLdSAcc, Z1);

        // 2) b = (I - A) mu
        // G_b = dL/db * d b / dA = -(dL/db) mu^T  (careful with shapes)
        DenseMatrix64F G_b = tmpK1; // k×k
        CommonOps.multTransB(-1.0, dLdb, mu, G_b);
        DenseMatrix64F Z2 = frechetHelper.frechetAdjointExp(blockDParams, G_b, t);
        CommonOps.scale(-t, Z2);
        CommonOps.addEquals(dLdSAcc, Z2);

        // 3) V = Sigma_stat - A Sigma_stat A^T
        // Contribution through A: G_V = -(dLdV A Sigma_stat + (dLdV)^T A Sigma_stat)
        DenseMatrix64F temp1 = tmpK2;
        DenseMatrix64F temp2 = tmpK3;
        DenseMatrix64F dLdV_T = tmpK4;

        CommonOps.mult(dLdV, A, temp1);
        CommonOps.mult(temp1, Sigma_stat, temp2);

        CommonOps.transpose(dLdV, dLdV_T);
        CommonOps.mult(dLdV_T, A, temp1);
        CommonOps.multAdd(temp1, Sigma_stat, temp2);

        DenseMatrix64F G_V = tmpK1;
        CommonOps.scale(-1.0, temp2, G_V);

        DenseMatrix64F Z3 = frechetHelper.frechetAdjointExp(blockDParams, G_V, t);
        CommonOps.scale(-t, Z3);
        CommonOps.addEquals(dLdSAcc, Z3);
    }

    // ---------- 3. Lyapunov contribution ----------

    private void computeLyapunovContribution(DenseMatrix64F Sigma_stat,
                                             DenseMatrix64F dLdSigmaStat,
                                             DenseMatrix64F dLdSAcc,
                                             double[] blockDParams) {
        lyapAdjHelper.accumulateLyapunovContribution(Sigma_stat, dLdSigmaStat, blockDParams, dLdSAcc);
    }
}
