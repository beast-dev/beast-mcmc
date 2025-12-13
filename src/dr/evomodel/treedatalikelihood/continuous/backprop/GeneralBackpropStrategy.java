package dr.evomodel.treedatalikelihood.continuous.backprop;

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import static dr.evomodel.treedatalikelihood.continuous.backprop.GeneralOUFrechetUtils.frechetAdjExp;

/**
 * General backpropagation strategy for OU models with a full,
 * unconstrained strength-of-selection matrix S.
 *
 * This implements the mathematically-correct Frechét adjoint derivatives
 * for all OU primitive parameters:
 *
 *   - dL/dS
 *   - dL/dSigma       (optional, derived via Lyapunov equation)
 *   - dL/dSigma_stat
 *   - dL/dMu
 *
 * It corresponds to the “GENERAL” case: S has no special structure
 * such as diagonal, block-diagonal, Schur form, or known eigenbasis.
 *
 * @author Filippo Monti
 */
public final class GeneralBackpropStrategy implements PrimitiveParameterBackpropStrategy {

    private final boolean debug;
    private final boolean computeSigmaGradient;

    // ---------------------------------------------------------------------
    // Constructors
    // ---------------------------------------------------------------------

    public GeneralBackpropStrategy() {
        this(false, false);
    }

    public GeneralBackpropStrategy(boolean debug, boolean computeSigmaGradient) {
        this.debug = debug;
        this.computeSigmaGradient = computeSigmaGradient;
    }

    // ---------------------------------------------------------------------
    // Main Backpropagation
    // ---------------------------------------------------------------------

    @Override
    public PrimitiveGradientSet backprop(
            ContinuousTraitBackpropGradient.OUBranchCache cache,
            MessageBackprop.Result leafGrads) {

        final int k = cache.A.numRows;

        if (debug) {
            System.err.println("\n=== GeneralBackpropStrategy Debug ===");
            printMatrix("dL/dA", leafGrads.dLdA);
            printMatrix("dL/dV", leafGrads.dLdV);
            printMatrix("dL/db", leafGrads.dLdb);
        }

        // ------------------------------------------------------------------
        // Allocate primitive gradients
        // ------------------------------------------------------------------

        DenseMatrix64F dLdS         = new DenseMatrix64F(k, k);
        DenseMatrix64F dLdSigmaStat = new DenseMatrix64F(k, k);
        DenseMatrix64F dLdMu        = new DenseMatrix64F(k, 1);

        // ------------------------------------------------------------------
        // 1. dL/dmu
        //    b = (I - A) mu  ⇒  dL/dμ = (I - A)^T · dL/db
        // ------------------------------------------------------------------

        DenseMatrix64F IminusA = new DenseMatrix64F(k, k);
        CommonOps.setIdentity(IminusA);
        CommonOps.subtractEquals(IminusA, cache.A);

        CommonOps.multTransA(IminusA, leafGrads.dLdb, dLdMu);

        // ------------------------------------------------------------------
        // 2. dL/dSigma_stat
        //    V = Σ_stat - A Σ_stat A^T
        //    ⇒ dL/dΣ_stat = dL/dV - A^T · dL/dV · A
        // ------------------------------------------------------------------

        DenseMatrix64F AT_dLdV     = new DenseMatrix64F(k, k);
        DenseMatrix64F AT_dLdV_A   = new DenseMatrix64F(k, k);

        CommonOps.multTransA(cache.A, leafGrads.dLdV, AT_dLdV);
        CommonOps.mult(AT_dLdV, cache.A, AT_dLdV_A);

        CommonOps.addEquals(dLdSigmaStat, leafGrads.dLdV);
        CommonOps.subtractEquals(dLdSigmaStat, AT_dLdV_A);

        if (debug) {
            printMatrix("dL/dμ", dLdMu);
            printMatrix("dL/dΣ_stat", dLdSigmaStat);
        }

        // ------------------------------------------------------------------
        // Precompute -S t
        // ------------------------------------------------------------------

        DenseMatrix64F minus_St = cache.S.copy();
        CommonOps.scale(-cache.t, minus_St);

        // ------------------------------------------------------------------
        // 3. Contribution 1: A = exp(-S t)
        // ------------------------------------------------------------------

        if (debug) System.err.println("\n--- Contribution 1 (A) ---");

        DenseMatrix64F Adj1 = new DenseMatrix64F(k, k);
        frechetAdjExp(minus_St, leafGrads.dLdA, Adj1);

        DenseMatrix64F contrib1 = Adj1.copy();
        CommonOps.scale(-cache.t, contrib1);
        CommonOps.addEquals(dLdS, contrib1);

        if (debug) {
            printMatrix("Adj1", Adj1);
            printMatrix("contrib1", contrib1);
        }

        // ------------------------------------------------------------------
        // 4. Contribution 2: b = (I-A) μ
        //    db/dA = -μ ⊗ 1^T
        // ------------------------------------------------------------------

        if (debug) System.err.println("\n--- Contribution 2 (b) ---");

        DenseMatrix64F muT = new DenseMatrix64F(1, k);
        CommonOps.transpose(cache.mu, muT);

        DenseMatrix64F G_b = new DenseMatrix64F(k, k);
        CommonOps.mult(leafGrads.dLdb, muT, G_b);
        CommonOps.scale(-1.0, G_b);

        DenseMatrix64F Adj2 = new DenseMatrix64F(k, k);
        frechetAdjExp(minus_St, G_b, Adj2);

        DenseMatrix64F contrib2 = Adj2.copy();
        CommonOps.scale(-cache.t, contrib2);
        CommonOps.addEquals(dLdS, contrib2);

        if (debug) {
            printMatrix("G_b", G_b);
            printMatrix("Adj2", Adj2);
            printMatrix("contrib2", contrib2);
        }

        // ------------------------------------------------------------------
        // 5. Contribution 3: V = Σ_stat - A Σ_stat A^T
        //
        // dV/dA[E] = -(E Σ_stat A^T + A Σ_stat E^T)
        //
        // G_V = -(dL/dV A Σ_stat + (dL/dV)^T A Σ_stat)
        // ------------------------------------------------------------------

        if (debug) System.err.println("\n--- Contribution 3 (V) ---");

        DenseMatrix64F dLdV_T = new DenseMatrix64F(k, k);
        CommonOps.transpose(leafGrads.dLdV, dLdV_T);

        DenseMatrix64F temp  = new DenseMatrix64F(k, k);
        DenseMatrix64F term1 = new DenseMatrix64F(k, k);
        DenseMatrix64F term2 = new DenseMatrix64F(k, k);

        CommonOps.mult(leafGrads.dLdV, cache.A, temp);
        CommonOps.mult(temp, cache.sigmaStat, term1);

        CommonOps.mult(dLdV_T, cache.A, temp);
        CommonOps.mult(temp, cache.sigmaStat, term2);

        DenseMatrix64F G_V = new DenseMatrix64F(k, k);
        CommonOps.add(term1, term2, G_V);
        CommonOps.scale(-1.0, G_V);

        DenseMatrix64F Adj3 = new DenseMatrix64F(k, k);
        frechetAdjExp(minus_St, G_V, Adj3);

        DenseMatrix64F contrib3 = Adj3.copy();
        CommonOps.scale(-cache.t, contrib3);
        CommonOps.addEquals(dLdS, contrib3);

        if (debug) {
            printMatrix("dL/dV^T", dLdV_T);
            printMatrix("term1", term1);
            printMatrix("term2", term2);
            printMatrix("G_V", G_V);
            printMatrix("Adj3", Adj3);
            printMatrix("contrib3", contrib3);
        }

        // ------------------------------------------------------------------
        // 6. Contribution 4: Σ_stat via Lyapunov equation
        //
        //     A Σ_stat + Σ_stat A^T = Σ
        //
        // dΣ_stat/dS[i,j]: computed elementwise via LyapunovSolver
        // ------------------------------------------------------------------

        if (debug) System.err.println("\n--- Contribution 4 (Σ_stat) ---");

        double[] S_flat = cache.S.getData();

        for (int i = 0; i < k; i++) {
            for (int j = 0; j < k; j++) {

                DenseMatrix64F dSigmaStat_dS_ij =
                        GeneralLyapunovSolver.derivativeWrtS(
                                S_flat, cache.sigmaStat, i, j, k
                        );

                double grad = frobeniusInner(dLdSigmaStat, dSigmaStat_dS_ij);

                dLdS.set(i, j, dLdS.get(i, j) + grad);

                if (debug && Math.abs(grad) > 1e-12) {
                    System.err.printf("  dS[%d,%d] contrib4 = %.6e%n", i, j, grad);
                }
            }
        }

        // ------------------------------------------------------------------
        // 7. Optional: convert dL/dΣ_stat → dL/dΣ via Lyapunov derivative
        // ------------------------------------------------------------------

        DenseMatrix64F dLdSigma = null;

        if (computeSigmaGradient) {
            if (debug) System.err.println("\n--- Computing dL/dΣ ---");

            dLdSigma = new DenseMatrix64F(k, k);

            for (int p = 0; p < k; p++) {
                for (int q = 0; q < k; q++) {

                    DenseMatrix64F dSigmaStat_dSigma_pq =
                            GeneralLyapunovSolver.derivativeWrtSigma(S_flat, p, q, k);

                    double grad = frobeniusInner(dLdSigmaStat, dSigmaStat_dSigma_pq);

                    dLdSigma.set(p, q, grad);
                }
            }

            // Σ is symmetric ⇒ symmetrize gradient
            for (int i = 0; i < k; i++) {
                for (int j = i + 1; j < k; j++) {
                    double s = 0.5 * (dLdSigma.get(i, j) + dLdSigma.get(j, i));
                    dLdSigma.set(i, j, s);
                    dLdSigma.set(j, i, s);
                }
            }

            if (debug) printMatrix("dL/dΣ", dLdSigma);
        }

        if (debug) {
            System.err.println("\n=== Final GeneralBackpropStrategy Gradients ===");
            printMatrix("dL/dS", dLdS);
            printMatrix("dL/dΣ_stat", dLdSigmaStat);
            printMatrix("dL/dμ", dLdMu);
        }

        return new PrimitiveGradientSet(dLdS, dLdSigma, dLdSigmaStat, dLdMu);
    }

    // ---------------------------------------------------------------------
    // Helpers
    // ---------------------------------------------------------------------

    private static double frobeniusInner(DenseMatrix64F A, DenseMatrix64F B) {
        double sum = 0.0;
        for (int i = 0; i < A.numRows; i++) {
            for (int j = 0; j < A.numCols; j++) {
                sum += A.get(i, j) * B.get(i, j);
            }
        }
        return sum;
    }

    private void printMatrix(String name, DenseMatrix64F m) {
        if (!debug) return;
        System.err.println(name + ":");
        for (int i = 0; i < m.numRows; i++) {
            System.err.print("  ");
            for (int j = 0; j < m.numCols; j++) {
                System.err.printf("%12.6f ", m.get(i, j));
            }
            System.err.println();
        }
    }
}
