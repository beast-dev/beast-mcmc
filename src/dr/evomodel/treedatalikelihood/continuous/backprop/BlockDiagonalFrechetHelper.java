package dr.evomodel.treedatalikelihood.continuous.backprop;

import dr.inference.model.BlockDiagonalCosSinMatrixParameter;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * Helper for Fréchet derivatives of exp(-St) when S = R D R^{-1},
 * with D block diagonal (1×1 and 2×2 blocks).
 *
 * Works entirely in the D-basis and uses BlockDiagonalFrechetIntegrator
 * for the 1×1–2×2 / 2×2–1×1 / 2×2–2×2 integrals.
 *
 * Valid for GENERAL invertible R (no orthogonality assumed).
 *
 * NOT thread-safe: owns internal workspaces.
 */

/**
 * @author Filippo Monti
 */
final class BlockDiagonalFrechetHelper {

    private final int k;
    private final int numBlocks;
    private final int[] blockStarts;
    private final int[] blockSizes;

    private final DenseMatrix64F R;      // change-of-basis
    private final DenseMatrix64F Rinv;   // R^{-1}

    // Dedicated D-basis workspaces (avoid aliasing with caller buffers)
    private final DenseMatrix64F dBasisTemp;
    private final DenseMatrix64F dBasisResult;

    // Workspace for compressed D^T params
    private final double[] blockDTransposeParams;

    private final double[] buf2a = new double[2];
    private final double[] buf4a = new double[4];
    private final double[] buf4b = new double[4];
    private final double[] buf4c = new double[4];

    BlockDiagonalFrechetHelper(BlockDiagonalCosSinMatrixParameter blockParam,
                               DenseMatrix64F R,
                               DenseMatrix64F Rinv) {

        this.k = blockParam.getRowDimension();
        this.blockStarts = blockParam.getBlockStarts();
        this.blockSizes  = blockParam.getBlockSizes();
        this.numBlocks   = blockParam.getNumBlocks();

        this.R    = R;
        this.Rinv = Rinv;

        this.dBasisTemp   = new DenseMatrix64F(k, k);
        this.dBasisResult = new DenseMatrix64F(k, k);

        this.blockDTransposeParams = new double[3 * k - 2];
    }

    /**
     * Compute L_exp^*(-St)[X] exploiting block diagonal D.
     *
     * Uses the identity:
     *
     *   L_exp^*(-St)[X]
     *     = R^{-T} · L_exp^*(-Dt)[ R^T X R^{-T} ] · R^T
     *
     * and the fact that in the D-basis
     *   L_exp^*(-Dt) = L_exp(-D^T t).
     */
    DenseMatrix64F frechetAdjointExp(double[] blockDParams,
                                     DenseMatrix64F X,
                                     double t) {

        final int diagLen = k;
        final int offLen  = k - 1;

        // Build compressed parameters for D^T
        System.arraycopy(blockDParams, 0,
                blockDTransposeParams, 0,
                diagLen);

        // upper(D^T) = lower(D)
        System.arraycopy(blockDParams, diagLen + offLen,
                blockDTransposeParams, diagLen,
                offLen);

        // lower(D^T) = upper(D)
        System.arraycopy(blockDParams, diagLen,
                blockDTransposeParams, diagLen + offLen,
                offLen);

        // (1) X̃ = R^T X R^{-T}
        DenseMatrix64F X_tilde = transformToDBasisAdjoint(X);

        // (2) forward Fréchet in D-basis with D^T
        DenseMatrix64F result_tilde =
                computeFrechetInDBasis(blockDTransposeParams, X_tilde, t);

        // (3) back: R^{-T} result̃ R^T
        return transformFromDBasisAdjoint(result_tilde);
    }

    /**
     * Compute L_exp(-St)[E] via D-basis and block-wise integrals.
     *
     * Uses:
     *   L_exp(-St)[E] = R · L_exp(-Dt)[ R^{-1} E R ] · R^{-1}
     */
    DenseMatrix64F frechetExp(double[] blockDParams,
                              DenseMatrix64F E,
                              double t) {

        // Step 1: Ẽ = R^{-1} E R
        DenseMatrix64F E_tilde = transformToDBasisForward(E);

        // Step 2: block-wise Fréchet in D-basis
        DenseMatrix64F result_tilde =
                computeFrechetInDBasis(blockDParams, E_tilde, t);

        // Step 3: back: R result_tilde R^{-1}
        return transformFromDBasisForward(result_tilde);
    }

    // ============================================================
    // Forward D-basis transforms
    // ============================================================

    private DenseMatrix64F transformToDBasisForward(DenseMatrix64F X) {
        // dBasisTemp = R^{-1} X
        CommonOps.mult(Rinv, X, dBasisTemp);
        // dBasisResult = (R^{-1} X) R
        CommonOps.mult(dBasisTemp, R, dBasisResult);
        return dBasisResult;
    }

    private DenseMatrix64F transformFromDBasisForward(DenseMatrix64F X_tilde) {
        // dBasisTemp = R X_tilde
        CommonOps.mult(R, X_tilde, dBasisTemp);
        // dBasisResult = (R X_tilde) R^{-1}
        CommonOps.mult(dBasisTemp, Rinv, dBasisResult);
        return dBasisResult;
    }

    // ============================================================
    // Adjoint D-basis transforms (NO explicit transposes allocated)
    // ============================================================

    private DenseMatrix64F transformToDBasisAdjoint(DenseMatrix64F X) {
        // dBasisTemp = R^T X
        CommonOps.multTransA(R, X, dBasisTemp);
        // dBasisResult = (R^T X) R^{-T}
        CommonOps.multTransB(dBasisTemp, Rinv, dBasisResult);
        return dBasisResult;
    }

    private DenseMatrix64F transformFromDBasisAdjoint(DenseMatrix64F X_tilde) {
        // dBasisTemp = R^{-T} X_tilde
        CommonOps.multTransA(Rinv, X_tilde, dBasisTemp);
        // dBasisResult = (R^{-T} X_tilde) R^T
        CommonOps.multTransB(dBasisTemp, R, dBasisResult);
        return dBasisResult;
    }

    // ============================================================
    // Block-wise Fréchet in D-basis
    // ============================================================

    private DenseMatrix64F computeFrechetInDBasis(double[] blockDParams,
                                                  DenseMatrix64F E_tilde,
                                                  double t) {
        DenseMatrix64F result = new DenseMatrix64F(k, k);

        for (int bi = 0; bi < numBlocks; bi++) {
            int startI = blockStarts[bi];
            int sizeI  = blockSizes[bi];

            for (int bj = 0; bj < numBlocks; bj++) {
                int startJ = blockStarts[bj];
                int sizeJ  = blockSizes[bj];

                computeBlockFrechetIntegral(blockDParams, E_tilde, t,
                        startI, sizeI, startJ, sizeJ, result);
            }
        }
        return result;
    }

    private void computeBlockFrechetIntegral(double[] blockDParams,
                                             DenseMatrix64F E_tilde,
                                             double t,
                                             int startI, int sizeI,
                                             int startJ, int sizeJ,
                                             DenseMatrix64F result) {

        final int upperOffset = k;
        final int lowerOffset = k + (k - 1);

        if (sizeI == 1 && sizeJ == 1) {
            double lambda_i = -t * blockDParams[startI];
            double lambda_j = -t * blockDParams[startJ];
            double e_ij     = E_tilde.get(startI, startJ);

            double integral;
            if (Math.abs(lambda_i - lambda_j) < 1e-10) {
                integral = e_ij * Math.exp(lambda_i);
            } else {
                integral = e_ij * (Math.exp(lambda_i) - Math.exp(lambda_j))
                        / (lambda_i - lambda_j);
            }
            result.set(startI, startJ, integral);
            return;
        }

        if (sizeI == 1 && sizeJ == 2) {
            double lambda_i = -t * blockDParams[startI];

            double a = -t * blockDParams[startJ];
            double d = -t * blockDParams[startJ + 1];
            double b = -t * blockDParams[upperOffset + startJ];
            double c = -t * blockDParams[lowerOffset + startJ];

            double e1 = E_tilde.get(startI, startJ);
            double e2 = E_tilde.get(startI, startJ + 1);

            BlockDiagonalFrechetIntegrator.integrate1x1_2x2(
                    lambda_i, a, b, c, d,
                    e1, e2,
                    buf2a,
                    buf4a
            );

            result.set(startI, startJ,     buf2a[0]);
            result.set(startI, startJ + 1, buf2a[1]);
            return;
        }

        if (sizeI == 2 && sizeJ == 1) {
            double a = -t * blockDParams[startI];
            double d = -t * blockDParams[startI + 1];
            double b = -t * blockDParams[upperOffset + startI];
            double c = -t * blockDParams[lowerOffset + startI];

            double lambda_j = -t * blockDParams[startJ];

            double e1 = E_tilde.get(startI,     startJ);
            double e2 = E_tilde.get(startI + 1, startJ);

            BlockDiagonalFrechetIntegrator.integrate2x2_1x1(
                    a, b, c, d,
                    lambda_j,
                    e1, e2,
                    buf2a,
                    buf4a
            );

            result.set(startI,     startJ, buf2a[0]);
            result.set(startI + 1, startJ, buf2a[1]);
            return;
        }

        if (sizeI == 2 && sizeJ == 2) {
            double ai = -t * blockDParams[startI];
            double di = -t * blockDParams[startI + 1];
            double bi = -t * blockDParams[upperOffset + startI];
            double ci = -t * blockDParams[lowerOffset + startI];

            double aj = -t * blockDParams[startJ];
            double dj = -t * blockDParams[startJ + 1];
            double bj = -t * blockDParams[upperOffset + startJ];
            double cj = -t * blockDParams[lowerOffset + startJ];

            double e11 = E_tilde.get(startI,     startJ);
            double e12 = E_tilde.get(startI,     startJ + 1);
            double e21 = E_tilde.get(startI + 1, startJ);
            double e22 = E_tilde.get(startI + 1, startJ + 1);

            BlockDiagonalFrechetIntegrator.integrate2x2_2x2(
                    ai, bi, ci, di,
                    aj, bj, cj, dj,
                    e11, e12, e21, e22,
                    buf4a,
                    buf4b,
                    buf4c
            );

            result.set(startI,     startJ,     buf4a[0]);
            result.set(startI,     startJ + 1, buf4a[1]);
            result.set(startI + 1, startJ,     buf4a[2]);
            result.set(startI + 1, startJ + 1, buf4a[3]);
        }
    }
}
