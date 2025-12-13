package dr.evomodel.treedatalikelihood.continuous.backprop;

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
* Backpropagation through the leaf message computation.
*
* Implements equations 28-37 from the backprop derivation:
*   - Canonical parameters (J, η, c) → branch quantities (A, V, b, r)
*
* Not thread-safe (owns internal workspaces). Create one instance per thread.
*/

/*
* @author Filippo Monti
*/

public final class MessageBackprop {

    private final boolean debug;

    // Cached dimension and reusable workspaces (reallocated only when k changes)
    private int kCached = -1;

    // k×k workspaces
    private DenseMatrix64F Wkk0;
    private DenseMatrix64F Wkk1;
    private DenseMatrix64F Wkk2;
    private DenseMatrix64F Wkk3;

    // k×1 workspaces
    private DenseMatrix64F Wk1_0;
    private DenseMatrix64F Wk1_1;

    public MessageBackprop(boolean debug) {
        this.debug = debug;
    }

    /** Result of backprop through leaf message. */
    public static final class Result {
        public final DenseMatrix64F dLdA;
        public final DenseMatrix64F dLdV;
        public final DenseMatrix64F dLdVInv;
        public final DenseMatrix64F dLdr;
        public final DenseMatrix64F dLdb;

        public Result(DenseMatrix64F dLdA, DenseMatrix64F dLdV, DenseMatrix64F dLdVInv,
                      DenseMatrix64F dLdr, DenseMatrix64F dLdb) {
            this.dLdA = dLdA;
            this.dLdV = dLdV;
            this.dLdVInv = dLdVInv;
            this.dLdr = dLdr;
            this.dLdb = dLdb;
        }
    }

    public Result backprop(ContinuousTraitBackpropGradient.OUBranchCache cache,
                           ContinuousTraitBackpropGradient.CanonicalAdjoint adjoint) {

        final DenseMatrix64F A = cache.A;
        final DenseMatrix64F V = cache.V;
        final DenseMatrix64F Vinv = cache.Vinv;
        final DenseMatrix64F r = cache.r;

        final DenseMatrix64F dLdJ = adjoint.dLdJ;
        final DenseMatrix64F dLdEta = adjoint.dLdEta;
        final double dLdC = adjoint.dLdC;

        final int k = A.numRows;
        ensureWorkspaces(k);

        // Allocate result matrices (kept as return objects; caller owns them)
        final DenseMatrix64F dLdA = new DenseMatrix64F(k, k);
        final DenseMatrix64F dLdV = new DenseMatrix64F(k, k);
        final DenseMatrix64F dLdVInv = new DenseMatrix64F(k, k);
        final DenseMatrix64F dLdr = new DenseMatrix64F(k, 1);
        final DenseMatrix64F dLdb = new DenseMatrix64F(k, 1);

        // Compute gradients
        computeGradientWrtA(A, Vinv, r, dLdJ, dLdEta, dLdA);
        computeGradientWrtVInv(A, V, r, dLdJ, dLdEta, dLdC, dLdVInv);
        symmetrizeInPlace(dLdVInv);
        computeGradientWrtV(Vinv, dLdVInv, dLdV);
        computeGradientWrtR(A, Vinv, r, dLdEta, dLdC, dLdr);
        computeGradientWrtB(dLdr, dLdb);

        return new Result(dLdA, dLdV, dLdVInv, dLdr, dLdb);
    }

    // -------------------------
    // Workspace management
    // -------------------------

    private void ensureWorkspaces(int k) {
        if (k == kCached) return;
        kCached = k;

        Wkk0 = new DenseMatrix64F(k, k);
        Wkk1 = new DenseMatrix64F(k, k);
        Wkk2 = new DenseMatrix64F(k, k);
        Wkk3 = new DenseMatrix64F(k, k);

        Wk1_0 = new DenseMatrix64F(k, 1);
        Wk1_1 = new DenseMatrix64F(k, 1);
    }

    private void zero(DenseMatrix64F M) {
        // CommonOps.fill exists in newer EJML; this is compatible with older DenseMatrix64F
        final double[] d = M.getData();
        for (int i = 0; i < d.length; i++) d[i] = 0.0;
    }

    // -------------------------
    // Gradients
    // -------------------------

    /**
     * Equation 28:
     *   ∂L/∂A = 2 V^{-1} A (∂L/∂J) + (V^{-1} r)(∂L/∂η)^T
     *
     * Avoids allocating v, VinvA, tmp, dLdEta^T, outer by reusing workspaces
     * and computing the outer product by loops.
     */
    private void computeGradientWrtA(DenseMatrix64F A, DenseMatrix64F Vinv, DenseMatrix64F r,
                                     DenseMatrix64F dLdJ, DenseMatrix64F dLdEta,
                                     DenseMatrix64F dLdA) {
        final int k = A.numRows;

        // v = Vinv * r   (k×1)
        final DenseMatrix64F v = Wk1_0;
        CommonOps.mult(Vinv, r, v);

        // Wkk0 = Vinv * A
        final DenseMatrix64F VinvA = Wkk0;
        CommonOps.mult(Vinv, A, VinvA);

        // dLdA = 2 * (VinvA * dLdJ)
        final DenseMatrix64F tmp = Wkk1;
        CommonOps.mult(VinvA, dLdJ, tmp);
        CommonOps.scale(2.0, tmp);
        dLdA.set(tmp);

        // dLdA += v * dLdEta^T  (outer product, O(k^2))
        for (int i = 0; i < k; i++) {
            final double vi = v.get(i, 0);
            for (int j = 0; j < k; j++) {
                dLdA.add(i, j, vi * dLdEta.get(j, 0));
            }
        }
    }

    /**
     * Equation 32:
     *   ∂L/∂V^{-1} =
     *       A (∂L/∂J) A^T
     *     + (A ∂L/∂η) r^T
     *     + (∂L/∂c) * (-1/2 rr^T + 1/2 V)
     *
     * Avoids:
     *  - allocating r^T
     *  - allocating term2 matrix (adds outer product by loops)
     *  - allocating V.copy() (adds 0.5*V by loops)
     */
    private void computeGradientWrtVInv(DenseMatrix64F A, DenseMatrix64F V, DenseMatrix64F r,
                                        DenseMatrix64F dLdJ, DenseMatrix64F dLdEta, double dLdC,
                                        DenseMatrix64F dLdVInv) {
        final int k = A.numRows;

        // dLdVInv = A dLdJ A^T
        final DenseMatrix64F AdLdJ = Wkk0;
        CommonOps.mult(A, dLdJ, AdLdJ);
        CommonOps.multTransB(AdLdJ, A, dLdVInv);

        // term2: (A dLdEta) r^T   -> add via outer product loops
        final DenseMatrix64F AdLdEta = Wk1_0;
        CommonOps.mult(A, dLdEta, AdLdEta);

        for (int i = 0; i < k; i++) {
            final double ai = AdLdEta.get(i, 0);
            for (int j = 0; j < k; j++) {
                dLdVInv.add(i, j, ai * r.get(j, 0));
            }
        }

        // term3: dLdC * (-1/2 rr^T + 1/2 V)
        // Build rrT into Wkk1: rrT = r r^T
        final DenseMatrix64F rrT = Wkk1;
        // rrT = r * r^T via loops (O(k^2)) to avoid r^T allocation
        for (int i = 0; i < k; i++) {
            final double ri = r.get(i, 0);
            for (int j = 0; j < k; j++) {
                rrT.set(i, j, ri * r.get(j, 0));
            }
        }

        // rrT = -0.5 * rrT
        CommonOps.scale(-0.5, rrT);

        // rrT += 0.5 * V  (no V.copy())
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < k; j++) {
                rrT.add(i, j, 0.5 * V.get(i, j));
            }
        }

        // dLdVInv += dLdC * rrT
        CommonOps.addEquals(dLdVInv, dLdC, rrT);
    }

    /** Equation 57: symmetrize matrix in place. */
    private void symmetrizeInPlace(DenseMatrix64F mat) {
        final int k = mat.numRows;
        for (int i = 0; i < k; i++) {
            for (int j = i + 1; j < k; j++) {
                final double sym = 0.5 * (mat.get(i, j) + mat.get(j, i));
                mat.set(i, j, sym);
                mat.set(j, i, sym);
            }
        }
    }

    /**
     * Equation 33:
     *   ∂L/∂V = -V^{-1} (∂L/∂V^{-1}) V^{-1}
     *
     * Reuses Wkk0 for tmp.
     */
    private void computeGradientWrtV(DenseMatrix64F Vinv, DenseMatrix64F dLdVInv,
                                     DenseMatrix64F dLdV) {
        final DenseMatrix64F tmp = Wkk0;
        CommonOps.mult(Vinv, dLdVInv, tmp);
        CommonOps.mult(tmp, Vinv, dLdV);
        CommonOps.scale(-1.0, dLdV);
    }

    /**
     * Equation 36:
     *   ∂L/∂r = V^{-T} A (∂L/∂η) - (∂L/∂c) V^{-1} r
     *
     * Uses multTransA to avoid forming V^{-T} explicitly.
     */
    private void computeGradientWrtR(DenseMatrix64F A, DenseMatrix64F Vinv, DenseMatrix64F r,
                                     DenseMatrix64F dLdEta, double dLdC,
                                     DenseMatrix64F dLdr) {
        final int k = A.numRows;

        // Wkk0 = Vinv^T * A
        final DenseMatrix64F VinvT_A = Wkk0;
        CommonOps.multTransA(Vinv, A, VinvT_A);

        // dLdr = (Vinv^T A) dLdEta
        CommonOps.mult(VinvT_A, dLdEta, dLdr);

        // v = Vinv * r
        final DenseMatrix64F v = Wk1_0;
        CommonOps.mult(Vinv, r, v);

        // dLdr += (-dLdC) * v
        CommonOps.addEquals(dLdr, -dLdC, v);

        if (debug) {
            // Optional lightweight sanity check (no allocations)
            if (dLdr.numRows != k || dLdr.numCols != 1) {
                throw new IllegalStateException("dLdr has wrong shape.");
            }
        }
    }

    /** Equation 37: ∂L/∂b = -∂L/∂r (no allocations). */
    private void computeGradientWrtB(DenseMatrix64F dLdr, DenseMatrix64F dLdb) {
        final int k = dLdr.numRows;
        for (int i = 0; i < k; i++) {
            dLdb.set(i, 0, -dLdr.get(i, 0));
        }
    }
}
