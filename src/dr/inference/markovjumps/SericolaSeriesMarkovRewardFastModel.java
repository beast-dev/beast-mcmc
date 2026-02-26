/*
 * SericolaSeriesMarkovRewardFastModel.java
 *
 * Reparametrized to accept BRANCH-NORMALIZED reward:
 *      Y = (total reward) / (branchLength)
 *
 * Relationship (densities):
 *   If f_R(x | t) is the density w.r.t. total reward x,
 *   then the density w.r.t. y is:
 *        f_Y(y | t) = t * f_R(y t | t)
 *
 */

package dr.inference.markovjumps;

import dr.evomodel.substmodel.DefaultEigenSystem;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.EigenSystem;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.Binomial;
import dr.math.GammaFunction;
import dr.math.matrixAlgebra.Vector;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Random;

public class SericolaSeriesMarkovRewardFastModel extends AbstractModel implements MarkovReward {

    private static final boolean DEBUG = false;

    // Dependencies
    private final SubstitutionModel underlyingSubstitutionModel;
    private final Parameter rewardRates;

    // Dimensions
    private final int dim;
    private final int dim2;
    private final int phi; // dim - 1

    // Accuracy / truncation
    private final double epsilon;
    private static final double X_EPS = 1e-15;

    // -----------------------
    // Sorting buffers / maps
    // -----------------------
    private final double[] unsortedQ;   // dim2
    private final double[] sortedQ;     // dim2
    private final double[] sortedR;     // dim

    // perm maps SORTED index -> ORIGINAL index
    private final int[] perm;
    // invPerm maps ORIGINAL index -> SORTED index
    private final int[] invPerm;

    // Mapping helpers for writing results directly to ORIGINAL order
    // outRowBaseBySorted[uS] = perm[uS] * dim
    private final int[] outRowBaseBySorted;
    // outColBySorted[vS] = perm[vS]
    private final int[] outColBySorted;

    private boolean permDirty = true;
    private boolean qDirty = true;
    private boolean cachesDirty = true;

    // -----------------------
    // Numerics (mutable)
    // -----------------------
    private double lambda;              // uniformization rate
    private final double[] P;           // dim2, P = I + sortedQ/lambda
    private final double[] invRewardDiff; // dim, invRewardDiff[h]=1/(sortedR[h]-sortedR[h-1]) for h=1..phi

    // Exponential for conditional probabilities (optional)
    private final EigenSystem eigenSystem;
    private EigenDecomposition eigenDecomposition;

    // Scratch for computeChnk Pn multiplication
    private final double[] PnScratch;   // dim2
    private final double[] matMulScratch; // dim2

    // ------------------------------------------------------------
    // Storage for C(h,n,k,uv) flattened into a single array
    //
    // Indexing:
    //   cOffset(h,n,k) = (((h * N1 + n) * N1 + k) * dim2)
    //   then add uv (in SORTED uv space)
    //
    // where N1 = capacityN + 1
    // ------------------------------------------------------------
    private double[] Cflat;

    // allocation capacity (does not imply validity)
    private int allocatedN = -1;
    private int allocatedN1 = 0;   // allocatedN + 1

    // computed extent valid for current numerics
    private int computedN = -1;    // -1 means nothing computed / invalid
    private double maxTime = 0.0;  // max time covered by current computed cache

    // Thread-safe scratch buffers to avoid allocations in computePdf hot path
    private static final class Scratch {
        int[] H;
        int[] NN;
        double[] lt;
        double[] premult;

        // precomputed per-t
        double[] xh;        // in [0,1]
        double[] oneMinus;  // 1-xh
        double[] ratio;     // xh/(1-xh) (only valid when interior)
        double[] w0;        // (1-xh)^n evolving across n
        boolean[] isZero;   // xh ~ 0
        boolean[] isOne;    // xh ~ 1

        double[] inc; // length dim2, in SORTED uv order
    }

    private final ThreadLocal<Scratch> tls = ThreadLocal.withInitial(Scratch::new);

    // Temporary objects for scalar wrapper (allocation-free)
    private final double[] scratchY1 = new double[1];
    private final double[] scratchT1 = new double[1];
    private final double[][] scratchW1 = new double[1][];

    // ---------------------------------------
    // MCMC store/restore bookkeeping
    // ---------------------------------------
    private boolean storedPermDirty, storedQDirty, storedCachesDirty;
    private boolean hasStoredState = false;

    // --------------------------------------------------------------------
    // Constructor (requested signature)
    // --------------------------------------------------------------------
    public SericolaSeriesMarkovRewardFastModel(
            SubstitutionModel underlyingSubstitutionModel,
            Parameter rewardRates,
            int dim,
            double epsilon
    ) {
        super("SericolaSeriesMarkovRewardFast");

        if (underlyingSubstitutionModel == null) {
            throw new IllegalArgumentException("underlyingSubstitutionModel must be non-null");
        }
        if (rewardRates == null) {
            throw new IllegalArgumentException("rewardRates must be non-null");
        }
        if (dim <= 0) {
            throw new IllegalArgumentException("dim must be > 0");
        }
        if (!(epsilon > 0.0 && epsilon < 1.0)) {
            throw new IllegalArgumentException("epsilon must be in (0,1)");
        }

        this.underlyingSubstitutionModel = underlyingSubstitutionModel;
        this.rewardRates = rewardRates;

        this.dim = dim;
        this.dim2 = dim * dim;
        this.phi = dim - 1;
        this.epsilon = epsilon;

        this.unsortedQ = new double[dim2];
        this.sortedQ = new double[dim2];
        this.sortedR = new double[dim];

        this.perm = new int[dim];
        this.invPerm = new int[dim];
        this.outRowBaseBySorted = new int[dim];
        this.outColBySorted = new int[dim];

        // Numerics buffers (allocated once, filled when dirty)
        this.P = new double[dim2];
        this.invRewardDiff = new double[dim];

        this.eigenSystem = new DefaultEigenSystem(dim);
        this.PnScratch = new double[dim2];
        this.matMulScratch = new double[dim2];

        // Register dependencies
        addModel(underlyingSubstitutionModel);
        addVariable(rewardRates);

        // Mark everything dirty initially
        this.permDirty = true;
        this.qDirty = true;
        this.cachesDirty = true;
    }

    // ============================================================
    // Public API: PDF (NOW EXPECTS y = reward/time)
    // ============================================================

    @Override
    public double[] computePdf(double y, double branchLength) {
        final double[] out = new double[dim2];
        computePdfIntoY(y, branchLength, out);
        return out;
    }

    @Override
    public double computePdf(double y, double branchLength, int i, int j) {
        throw new UnsupportedOperationException("Not implemented yet.");
    }

    public void computePdfIntoY(double y, double branchLength, double[] out) {
        if (out == null || out.length != dim2) {
            throw new IllegalArgumentException("out must be length dim*dim=" + dim2);
        }
        scratchY1[0] = y;
        scratchT1[0] = branchLength;
        scratchW1[0] = out;
        computePdfIntoY(scratchY1, scratchT1, false, scratchW1);
    }

    public void computePdfIntoY(double[] Y, double time, double[][] W) {
        computePdfIntoY(Y, new double[]{time}, false, W);
    }

    public void computePdfIntoY(double[] Y, double[] times, double[][] W) {
        computePdfIntoY(Y, times, false, W);
    }

    public void computePdfIntoY(double[] Y, double[] times, boolean parsimonious, double[][] W) {
        validateInputsY(Y, times, W);
        final int T = Y.length;
        if (T == 0) return;

        ensureNumericsUpToDate();

        validateAndZeroOutputs(W, T);

        final Scratch s = tls.get();
        ensureScratchCapacity(s, T);

        final double maxT = precomputePdfScratchY(Y, times, parsimonious, s);

        // Ensure C is available up to required N (pdf uses n+1)
        ensureCForTime(maxT, /*extraN=*/1);
        final int N = computedN - 1;

        accumulatePdfOverN_Y(W, s, T, N);
    }

    // ----------------------------------------------------------------
    // Backward-compat wrappers (if any code still passes total rewards)
    // These convert X -> Y = X/t and call computePdfIntoY.
    // ----------------------------------------------------------------
    @Deprecated
    public void computePdfInto(double reward, double branchLength, double[] out) {
        if (branchLength <= 0.0) throw new IllegalArgumentException("branchLength must be > 0");
        computePdfIntoY(reward / branchLength, branchLength, out);
    }

    @Deprecated
    public void computePdfInto(double[] X, double time, double[][] W) {
        if (time <= 0.0) throw new IllegalArgumentException("time must be > 0");
        final double[] Y = new double[X.length];
        for (int i = 0; i < X.length; i++) Y[i] = X[i] / time;
        computePdfIntoY(Y, time, W);
    }

    @Deprecated
    public void computePdfInto(double[] X, double[] times, double[][] W) {
        final double[] Y = new double[X.length];
        if (times.length == 1) {
            final double t = times[0];
            if (t <= 0.0) throw new IllegalArgumentException("time must be > 0");
            for (int i = 0; i < X.length; i++) Y[i] = X[i] / t;
        } else {
            if (times.length != X.length) throw new IllegalArgumentException("times length mismatch");
            for (int i = 0; i < X.length; i++) {
                final double t = times[i];
                if (t <= 0.0) throw new IllegalArgumentException("time must be > 0 at index " + i);
                Y[i] = X[i] / t;
            }
        }
        computePdfIntoY(Y, times, false, W);
    }

    // ============================================================
    // Input validation / utilities
    // ============================================================

    private void validateInputsY(double[] Y, double[] times, double[][] W) {
        if (Y == null || times == null || W == null) {
            throw new IllegalArgumentException("Y/times/W must be non-null");
        }
        final int T = Y.length;

        final boolean singleTime = (times.length == 1);
        if (!singleTime && times.length != T) {
            throw new IllegalArgumentException("Either times.length==1 or times.length==Y.length");
        }
        if (W.length != T) {
            throw new IllegalArgumentException("W.length must equal Y.length");
        }
    }

    private void validateAndZeroOutputs(double[][] W, int T) {
        for (int t = 0; t < T; ++t) {
            final double[] row = W[t];
            if (row == null || row.length != dim2) {
                throw new IllegalArgumentException("W[" + t + "] must be non-null and length dim*dim=" + dim2);
            }
            Arrays.fill(row, 0.0);
        }
    }

    private double precomputePdfScratchY(double[] Y, double[] times, boolean parsimonious, Scratch s) {
        final int T = Y.length;
        final boolean singleTime = (times.length == 1);

        double maxT = 0.0;

        if (singleTime) {
            final double time = times[0];
            if (time <= 0.0) throw new IllegalArgumentException("time must be > 0");
            maxT = time;

            final double lt = lambda * time;
            final double premult0 = Math.exp(-lt);
            final int NN0 = (parsimonious ? determineNumberOfSteps(time) : Integer.MAX_VALUE);

            for (int t = 0; t < T; ++t) {
                s.lt[t] = lt;
                s.premult[t] = premult0;
                s.NN[t] = NN0;

                final int h = getHfromY(Y[t]);
                s.H[t] = h;

                fillXhPrecompY(s, t, Y[t], h);
            }
            return maxT;
        }

        // per-time
        for (int t = 0; t < T; ++t) {
            final double time = times[t];
            if (time <= 0.0) throw new IllegalArgumentException("time must be > 0 at index " + t);
            if (time > maxT) maxT = time;

            final double lt = lambda * time;
            s.lt[t] = lt;
            s.premult[t] = Math.exp(-lt);
            s.NN[t] = (parsimonious ? determineNumberOfSteps(time) : Integer.MAX_VALUE);

            final int h = getHfromY(Y[t]);
            s.H[t] = h;

            fillXhPrecompY(s, t, Y[t], h);
        }
        return maxT;
    }

    private void fillXhPrecompY(Scratch s, int t, double y, int h) {
        final double invDiff = invRewardDiff[h];

        // xh = (y - r_{h-1}) / (r_h - r_{h-1})
        double xh = (y - sortedR[h - 1]) * invDiff;

        if (xh <= 0.0) {
            s.xh[t] = 0.0;
            s.isZero[t] = true;
            s.isOne[t]  = false;
            return;
        }
        if (xh >= 1.0) {
            s.xh[t] = 1.0;
            s.isZero[t] = false;
            s.isOne[t]  = true;
            return;
        }

        s.xh[t] = xh;
        s.isZero[t] = false;
        s.isOne[t]  = false;

        final double om = 1.0 - xh;
        s.oneMinus[t] = om;
        s.ratio[t]    = xh / om;
        s.w0[t]       = 1.0;
    }

    private void accumulatePdfOverN_Y(double[][] W, Scratch s, int T, int N) {
        for (int n = 0; n <= N; ++n) {

            for (int t = 0; t < T; ++t) {
                if (s.NN[t] < n) continue;

                final int h = s.H[t];

                // Jacobian factor: f_Y(y|t) = t f_R(yt|t), and here t = lt/lambda
                final double time = s.lt[t] / lambda;
                final double base = (lambda * invRewardDiff[h]) * s.premult[t] * time;

                if (s.isZero[t]) {
                    addDiffBlockToOriginalOrder(W[t], base,
                            cOffset(h, n + 1, 1),
                            cOffset(h, n + 1, 0));
                } else if (s.isOne[t]) {
                    addDiffBlockToOriginalOrder(W[t], base,
                            cOffset(h, n + 1, n + 1),
                            cOffset(h, n + 1, n));
                } else {
                    loopCyclePdfAddToOriginalOrderFast_Y(
                            h, n,
                            s.ratio[t],
                            s.w0[t],
                            s.premult[t],
                            time,
                            W[t],
                            s.inc
                    );
                }
            }

            // premult_{n+1}(t) = premult_n(t) * (lambda*t)/(n+1) = premult_n(t) * lt/(n+1)
            final double inv = 1.0 / (n + 1.0);
            for (int t = 0; t < T; ++t) {
                s.premult[t] *= s.lt[t] * inv;
            }

            // w0: (1-xh)^n -> (1-xh)^(n+1) for interior points
            for (int t = 0; t < T; ++t) {
                if (!s.isZero[t] && !s.isOne[t]) {
                    s.w0[t] *= s.oneMinus[t];
                }
            }
        }
    }

    private void loopCyclePdfAddToOriginalOrderFast_Y(
            int h, int n,
            double ratio,
            double w0,                 // (1-xh)^n
            double premult,
            double time,               // Jacobian
            double[] WtOriginal,
            double[] incSorted) {

        final double temp = (lambda * invRewardDiff[h]) * premult * time;

        Arrays.fill(incSorted, 0.0);

        // k = 0
        {
            final int aOff = cOffset(h, n + 1, 1);
            final int bOff = cOffset(h, n + 1, 0);
            for (int uv = 0; uv < dim2; ++uv) {
                incSorted[uv] += w0 * (Cflat[aOff + uv] - Cflat[bOff + uv]);
            }
        }

        // k = 1..n via Bernstein recurrence
        double w = w0;
        for (int k = 0; k < n; ++k) {
            w *= ((double) (n - k) / (double) (k + 1)) * ratio;
            final int kp1 = k + 1;

            final int aOff = cOffset(h, n + 1, kp1 + 1);
            final int bOff = cOffset(h, n + 1, kp1);

            for (int uv = 0; uv < dim2; ++uv) {
                incSorted[uv] += w * (Cflat[aOff + uv] - Cflat[bOff + uv]);
            }
        }

        // write to ORIGINAL order
        for (int uS = 0; uS < dim; ++uS) {
            final int outRowBase = outRowBaseBySorted[uS];
            final int inRowBase = uS * dim;
            for (int vS = 0; vS < dim; ++vS) {
                WtOriginal[outRowBase + outColBySorted[vS]] += temp * incSorted[inRowBase + vS];
            }
        }
    }

    // ============================================================
    // Public API: CDF (NOW EXPECTS y)
    // ============================================================

    @Override
    public double computeCdf(double y, double time, int i, int j) {
        return computeCdfY(new double[]{y}, time)[0][i * dim + j];
    }

    public double[] computeCdf(double y, double time) {
        return computeCdfY(new double[]{y}, time)[0];
    }

    public double[][] computeCdfY(double[] Y, double time) {
        if (time <= 0.0) throw new IllegalArgumentException("time must be > 0");

        ensureNumericsUpToDate();
        ensureCForTime(time, /*extraN=*/0);
        final int N = computedN;

        final int T = Y.length;
        final double[][] W = new double[T][dim2];
        final int[] H = new int[T];
        for (int t = 0; t < T; t++) H[t] = getHfromY(Y[t]);

        for (int n = 0; n <= N; ++n) {
            accumulateCdfAddToOriginalOrder_Y(W, Y, H, n, time);
        }
        return W;
    }

    // Backward-compat CDF wrappers: X -> Y
    @Deprecated
    public double[][] computeCdf(double[] X, double time) {
        if (time <= 0.0) throw new IllegalArgumentException("time must be > 0");
        final double[] Y = new double[X.length];
        for (int i = 0; i < X.length; i++) Y[i] = X[i] / time;
        return computeCdfY(Y, time);
    }

    // ============================================================
    // Inner kernels
    // ============================================================

    /**
     * Hot inner kernel for PDF: add contribution for given block directly into Wt (ORIGINAL order).
     *
     * The cache Cflat is indexed by SORTED uv, but we write into ORIGINAL uv via precomputed maps:
     *   sorted (uS,vS) -> original index (perm[uS], perm[vS]).
     */
    private void addDiffBlockToOriginalOrder(double[] WtOriginal, double temp, int aOff, int bOff) {
        // Adds temp*(C[a]-C[b]) for all uv, mapping sorted uv -> original uv
        for (int uS = 0; uS < dim; ++uS) {
            final int outRowBase = outRowBaseBySorted[uS];
            final int inRowBase = uS * dim;
            for (int vS = 0; vS < dim; ++vS) {
                final int uvS = inRowBase + vS;
                final int outIdx = outRowBase + outColBySorted[vS];
                WtOriginal[outIdx] += temp * (Cflat[aOff + uvS] - Cflat[bOff + uvS]);
            }
        }
    }

    // CDF accumulation (writes in ORIGINAL order)
    private void accumulateCdfAddToOriginalOrder_Y(double[][] W, double[] Y, int[] H, int n, double time) {

        final double lt = lambda * time;
        final double premult = Math.exp(
                -lt + n * (Math.log(lambda) + Math.log(time)) - GammaFunction.lnGamma(n + 1.0)
        );

        for (int t = 0; t < Y.length; ++t) {
            final double y = Y[t];
            final int h = H[t];

            double xh = (y - sortedR[h - 1]) / (sortedR[h] - sortedR[h - 1]);

            // Guard numerical drift
            if (xh <= 0.0) xh = 0.0;
            else if (xh >= 1.0) xh = 1.0;

            // Compute increment in SORTED uv space, then write to ORIGINAL
            final double[] incSorted = new double[dim2];
            for (int k = 0; k <= n; k++) {
                final double binomialCoef = Binomial.choose(n, k) * Math.pow(xh, k) * Math.pow(1.0 - xh, n - k);
                final int cOff = cOffset(h, n, k);
                for (int uv = 0; uv < dim2; ++uv) {
                    incSorted[uv] += binomialCoef * Cflat[cOff + uv];
                }
            }

            final double[] WtOriginal = W[t];
            for (int uS = 0; uS < dim; ++uS) {
                final int outRowBase = outRowBaseBySorted[uS];
                final int inRowBase = uS * dim;
                for (int vS = 0; vS < dim; ++vS) {
                    final int outIdx = outRowBase + outColBySorted[vS];
                    WtOriginal[outIdx] += premult * incSorted[inRowBase + vS];
                }
            }
        }
    }

    // ============================================================
    // C cache growth + computation
    // ============================================================

    private void ensureCForTime(double time, int extraN) {
        if (time <= 0.0) throw new IllegalArgumentException("time must be > 0");

        final int requiredN = determineNumberOfSteps(time) + extraN;

        // ensure allocation first
        ensureCapacityN(requiredN);

        // track max time covered by computed cache
        if (time > maxTime) maxTime = time;

        // compute (or recompute) if current cache does not cover requiredN
        if (computedN < requiredN) {
            computeChnk(requiredN);     // compute up to requiredN
            computedN = requiredN;
        }
    }

    private void ensureCapacityN(int requiredN) {
        int newAlloc = (allocatedN < 0) ? requiredN : allocatedN;
        while (newAlloc < requiredN) {
            newAlloc = Math.max(requiredN, (int) (newAlloc * 1.5) + 1);
        }

        if (newAlloc <= allocatedN && Cflat != null) return;

        final int newN1 = newAlloc + 1;
        final long totalDoubles = (long) (phi + 1) * (long) newN1 * (long) newN1 * (long) dim2;

        if (totalDoubles > Integer.MAX_VALUE) {
            throw new IllegalStateException("C array too large: " + totalDoubles + " doubles");
        }

        Cflat = new double[(int) totalDoubles];
        allocatedN = newAlloc;
        allocatedN1 = newN1;

        // Allocation changed => any previously computed content is meaningless under new layout
        computedN = -1;
        maxTime = 0.0;
    }

    /**
     * Compute all C(h,n,k,uv) up to N using Sericola recursions,
     * assuming P, lambda, sortedR are current and consistent.
     *
     * IMPORTANT: Works in SORTED uv space.
     */
    private void computeChnk(int N) {
        if (Cflat == null || allocatedN < 0) {
            throw new IllegalStateException("C storage not initialized");
        }
        clearCflatUpToN(N);

        // Pn = I
        Arrays.fill(PnScratch, 0.0);
        for (int u = 0; u < dim; ++u) {
            PnScratch[idx(u, u)] = 1.0;
        }

        // Corner cases
        for (int h = 1; h <= phi; ++h) {
            final int off = cOffset(h, 0, 0);
            for (int u = 0; u <= h - 1; ++u) {
                Cflat[off + idx(u, u)] = 1.0;
            }
        }

        for (int n = 1; n <= N; ++n) {

            // Forward sweep
            for (int h = 1; h <= phi; ++h) {
                for (int k = 1; k <= n; ++k) {

                    for (int u = h; u <= phi; ++u) {
                        final int uvRow = u * dim;

                        final double cScalar = (sortedR[u] - sortedR[h]) / (sortedR[u] - sortedR[h - 1]);
                        final double dScalar = (sortedR[h] - sortedR[h - 1]) / (sortedR[u] - sortedR[h - 1]);

                        final int curOff = cOffset(h, n, k);
                        final int curOffKm1 = cOffset(h, n, k - 1);
                        final int prevOff = cOffset(h, n - 1, k - 1);

                        for (int v = 0; v <= phi; ++v) {
                            double cVal = cScalar * Cflat[curOffKm1 + (uvRow + v)];

                            double dVal = 0.0;
                            final int pRow = uvRow;
                            for (int w = 0; w <= phi; ++w) {
                                dVal += P[pRow + w] * Cflat[prevOff + (w * dim + v)];
                            }
                            dVal *= dScalar;

                            Cflat[curOff + (uvRow + v)] = cVal + dVal;
                        }
                    }
                }

                if (h < phi) {
                    final int srcOff = cOffset(h, n, n);
                    final int dstOff = cOffset(h + 1, n, 0);
                    for (int u = h + 1; u <= phi; ++u) {
                        System.arraycopy(Cflat, srcOff + u * dim, Cflat, dstOff + u * dim, dim);
                    }
                }
            }

            // Pn = Pn * P
            rightMultiply(PnScratch, P, matMulScratch);
            System.arraycopy(matMulScratch, 0, PnScratch, 0, dim2);

            // C(phi,n,n) = Pn for u=0..phi-1
            {
                final int off = cOffset(phi, n, n);
                for (int u = 0; u <= phi - 1; ++u) {
                    System.arraycopy(PnScratch, u * dim, Cflat, off + u * dim, dim);
                }
            }

            // Backward sweep
            for (int h = phi; h >= 1; --h) {
                for (int k = n - 1; k >= 0; --k) {

                    for (int u = 0; u <= h - 1; u++) {
                        final int uvRow = u * dim;

                        final double cScalar = (sortedR[h - 1] - sortedR[u]) / (sortedR[h] - sortedR[u]);
                        final double dScalar = (sortedR[h] - sortedR[h - 1]) / (sortedR[h] - sortedR[u]);

                        final int curOff = cOffset(h, n, k);
                        final int curOffKp1 = cOffset(h, n, k + 1);
                        final int prevOff = cOffset(h, n - 1, k);

                        for (int v = 0; v <= phi; ++v) {
                            double cVal = cScalar * Cflat[curOffKp1 + (uvRow + v)];

                            double dVal = 0.0;
                            final int pRow = uvRow;
                            for (int w = 0; w <= phi; ++w) {
                                dVal += P[pRow + w] * Cflat[prevOff + (w * dim + v)];
                            }
                            dVal *= dScalar;

                            Cflat[curOff + (uvRow + v)] = cVal + dVal;
                        }
                    }

                    if (h >= 2) {
                        final int srcOff = cOffset(h, n, 0);
                        final int dstOff = cOffset(h - 1, n, n);
                        for (int u = 0; u <= h - 2; ++u) {
                            System.arraycopy(Cflat, srcOff + u * dim, Cflat, dstOff + u * dim, dim);
                        }
                    }
                }
            }
        }
    }

    private void rightMultiply(double[] A, double[] B, double[] out) {
        Arrays.fill(out, 0.0);
        for (int i = 0; i < dim; ++i) {
            final int ioff = i * dim;
            for (int k = 0; k < dim; ++k) {
                final double aik = A[ioff + k];
                final int koff = k * dim;
                for (int j = 0; j < dim; ++j) {
                    out[ioff + j] += aik * B[koff + j];
                }
            }
        }
    }

    // ============================================================
    // Dirty-handling: permutation, sortedQ, numerics
    // ============================================================

    private void ensureNumericsUpToDate() {
        ensurePermutationUpToDate();
        ensureSortedQUpToDate();

        if (!cachesDirty) return;

        lambda = determineLambda(sortedQ);
        if (!(lambda > 0.0) || Double.isNaN(lambda) || Double.isInfinite(lambda)) {
            throw new IllegalStateException("Invalid uniformization rate lambda=" + lambda);
        }
        fillP(sortedQ, lambda, P);

        Arrays.fill(invRewardDiff, 0.0);
        for (int h = 1; h < dim; h++) {
            final double d = sortedR[h] - sortedR[h - 1];
            if (!(d > 0.0)) {
                throw new IllegalArgumentException("rewardRates must be strictly increasing after sorting; tie at h=" + h);
            }
            invRewardDiff[h] = 1.0 / d;
        }

        eigenDecomposition = null;

        invalidateCComputedExtent(); // this sets computedN=-1 and maxTime=0
        cachesDirty = false;
    }

    private void ensurePermutationUpToDate() {
        if (!permDirty) return;

        final double[] rewardVals = rewardRates.getParameterValues();
        if (rewardVals.length != dim) {
            throw new IllegalArgumentException("rewardRates length mismatch: expected " + dim + " got " + rewardVals.length);
        }

        Integer[] idx = new Integer[dim];
        for (int i = 0; i < dim; i++) idx[i] = i;
        Arrays.sort(idx, Comparator.comparingDouble(i -> rewardVals[i]));

        for (int s = 0; s < dim; s++) {
            perm[s] = idx[s];            // sorted index -> original index
            sortedR[s] = rewardVals[perm[s]];
        }
        for (int o = 0; o < dim; o++) {
            invPerm[o] = -1;
        }
        for (int s = 0; s < dim; s++) {
            invPerm[perm[s]] = s;        // original -> sorted
        }

        // Build write-back maps for sorted (uS,vS) -> original index
        for (int uS = 0; uS < dim; uS++) {
            outRowBaseBySorted[uS] = perm[uS] * dim;
            outColBySorted[uS] = perm[uS];
        }

        permDirty = false;
        qDirty = true;
        cachesDirty = true;

        if (DEBUG) {
            System.err.println("[DEBUG] rewardRates (original): " +
                    Arrays.toString(rewardRates.getParameterValues()));
            System.err.println("[DEBUG] sortedR: " + Arrays.toString(sortedR));
            System.err.println("[DEBUG] perm (sorted -> original): " + Arrays.toString(perm));
            System.err.println("[DEBUG] invPerm (original -> sorted): " + Arrays.toString(invPerm));
        }
    }

    private void ensureSortedQUpToDate() {
        if (!qDirty) return;

        underlyingSubstitutionModel.getInfinitesimalMatrix(unsortedQ);
        if (unsortedQ.length != dim2) {
            throw new IllegalStateException("Internal unsortedQ length mismatch");
        }

        for (int iS = 0; iS < dim; ++iS) {
            final int iO = perm[iS];          // original row
            final int rowS = iS * dim;
            final int rowO = iO * dim;
            for (int jS = 0; jS < dim; ++jS) {
                final int jO = perm[jS];      // original col
                sortedQ[rowS + jS] = unsortedQ[rowO + jO];
            }
        }

        qDirty = false;
        cachesDirty = true;
    }

    // ============================================================
    // Helpers
    // ============================================================

    private void ensureScratchCapacity(Scratch s, int T) {
        if (s.H == null || s.H.length < T) s.H = new int[T];
        if (s.NN == null || s.NN.length < T) s.NN = new int[T];
        if (s.lt == null || s.lt.length < T) s.lt = new double[T];
        if (s.premult == null || s.premult.length < T) s.premult = new double[T];

        if (s.xh == null || s.xh.length < T) s.xh = new double[T];
        if (s.oneMinus == null || s.oneMinus.length < T) s.oneMinus = new double[T];
        if (s.ratio == null || s.ratio.length < T) s.ratio = new double[T];
        if (s.w0 == null || s.w0.length < T) s.w0 = new double[T];
        if (s.isZero == null || s.isZero.length < T) s.isZero = new boolean[T];
        if (s.isOne == null || s.isOne.length < T) s.isOne = new boolean[T];

        if (s.inc == null || s.inc.length != dim2) s.inc = new double[dim2];
    }

    private int idx(int i, int j) {
        return i * dim + j;
    }

    private int cOffset(int h, int n, int k) {
        return (((h * allocatedN1 + n) * allocatedN1 + k) * dim2);
    }

    private void clearCflatUpToN(int N) {
        // We will write/read C(h,n,k,uv) for h=0..phi and n,k=0..N.
        // Clear exactly the range that can be touched.
        final int end = cOffset(phi, N, N) + dim2; // exclusive end index
        Arrays.fill(Cflat, 0, end, 0.0);
    }

    private double determineLambda(double[] QrowMajor) {
        double minDiag = QrowMajor[0];          // Q[0,0]
        for (int i = 1; i < dim; ++i) {
            final double d = QrowMajor[i * dim + i];
            if (d < minDiag) minDiag = d;
        }
        return -minDiag;
    }

    private void invalidateCComputedExtent() {
        computedN = -1;
        maxTime = 0.0;
        // do NOT touch allocation (Cflat, allocatedN, allocatedN1)
    }

    private void fillP(double[] Qsorted, double lambda, double[] Pout) {
        final double inv = 1.0 / lambda;
        for (int i = 0; i < dim; ++i) {
            final int ioff = i * dim;
            for (int j = 0; j < dim; ++j) {
                Pout[ioff + j] = (i == j ? 1.0 : 0.0) + Qsorted[ioff + j] * inv;
            }
        }
    }

    final Random rng = new Random(666);

    /**
     * Determine h based on y in SORTED reward space.
     * Returns an h in [1..phi].
     */
    private int getHfromY(double y) {

        final double lo = sortedR[0];
        final double hi = sortedR[phi];

        if (y < lo - 1e-10) {
            throw new IllegalArgumentException("y must be >= minRewardRate");
        }

        if (y > hi + 1e-10) {
            throw new IllegalArgumentException("y must be <= maxRewardRate");
        }

        int idx = Arrays.binarySearch(sortedR, 0, phi + 1, y);

        // Exact knot
        if (idx >= 0) {
            if (idx == 0) return 1;        // only possible interval
            if (idx == phi) return phi;    // only possible interval
            return idx + 1;                // deterministic: upper interval
        }

        // Strictly inside an interval
        int insertionPoint = -idx - 1;   // in [1..phi]
        return Math.max(1, Math.min(phi, insertionPoint));
    }

    /**
     * Fast Poisson CDF truncation:
     * find smallest i such that sum_{k=0}^i e^{-lt} (lt)^k / k! >= 1-epsilon.
     */
    private int determineNumberOfSteps(double time) {
        final double target = 1.0 - epsilon;
        final double lt = lambda * time;

        double term = Math.exp(-lt);
        double sum = term;

        int i = 0;
        final int hardCap = Math.max(5000, (int) (lt + 10.0 * Math.sqrt(lt + 1.0)));

        while (sum < target && i < hardCap) {
            i++;
            term *= lt / i;
            sum += term;
        }
        return i;
    }

    private double[][] squareMatrix(final double[] matRowMajor) {
        final double[][] rtn = new double[dim][dim];
        for (int i = 0; i < dim; ++i) {
            System.arraycopy(matRowMajor, i * dim, rtn[i], 0, dim);
        }
        return rtn;
    }

    // ============================================================
    // Conditional probabilities (matrix exponential of sortedQ)
    // ============================================================

    private EigenDecomposition getEigenDecomposition() {
        ensureNumericsUpToDate();
        if (eigenDecomposition == null) {
            eigenDecomposition = eigenSystem.decomposeMatrix(squareMatrix(sortedQ));
        }
        return eigenDecomposition;
    }

    public double[] computeConditionalProbabilities(double distance) {
        ensureNumericsUpToDate();
        final double[] matrix = new double[dim2];
        eigenSystem.computeExponential(getEigenDecomposition(), distance, matrix);
        return matrix;
    }

    public double computeConditionalProbability(double distance, int i, int j) {
        ensureNumericsUpToDate();
        return eigenSystem.computeExponential(getEigenDecomposition(), distance, i, j);
    }

    // ============================================================
    // BEAST Model callbacks
    // ============================================================

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == underlyingSubstitutionModel) {
            qDirty = true;
            cachesDirty = true;
            fireModelChanged();
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == rewardRates) {
            permDirty = true;
            qDirty = true;
            cachesDirty = true;
            fireModelChanged();
        }
    }

    @Override
    protected void storeState() {
        storedPermDirty = permDirty;
        storedQDirty = qDirty;
        storedCachesDirty = cachesDirty;
        hasStoredState = true;
    }

    @Override
    protected void restoreState() {
        if (!hasStoredState) return;

        permDirty = true;
        qDirty = true;
        cachesDirty = true;

        invalidateCComputedExtent();
        eigenDecomposition = null;
    }

    @Override
    protected void acceptState() {
        // nothing
    }

    // ============================================================
    // Debug/info
    // ============================================================

    @Override
    public String toString() {
        ensureNumericsUpToDate();
        StringBuilder sb = new StringBuilder();
        sb.append("sortedQ: ").append(new Vector(sortedQ)).append("\n");
        sb.append("sortedR: ").append(new Vector(sortedR)).append("\n");
        sb.append("lambda: ").append(lambda).append("\n");
        sb.append("allocateN: ").append(allocatedN).append("\n");
        sb.append("maxTime: ").append(maxTime).append("\n");
        sb.append("cprob at maxTime: ").append(new Vector(computeConditionalProbabilities(maxTime))).append("\n");
        return sb.toString();
    }

    // ============================================================
    // NO-ALLOCATION derivative of reward-density w.r.t. NORMALIZED REWARD y
    // ============================================================
    //
    // Computes: d/dy f_{ij}(y | t)  where y = totalReward / t
    //
    // Convention for boundary snapping (xh==0 or 1):
    //   returns 0 derivative (consistent with your "snap-to-boundary" semantics).

    public void computePdfDerivativeWrtYInto(double y, double branchLength, double[] out) {
        if (out == null || out.length != dim2) {
            throw new IllegalArgumentException("out must be length dim*dim=" + dim2);
        }
        if (branchLength <= 0.0) {
            throw new IllegalArgumentException("branchLength must be > 0");
        }
        Arrays.fill(out, 0.0);

        ensureNumericsUpToDate();

        // Determine h from y
        final int h = getHfromY(y);

        // Precompute xh and boundary flags
        final Scratch s = tls.get();
        ensureScratchCapacity(s, 1);
        fillXhPrecompY(s, 0, y, h);

        if (s.isZero[0] || s.isOne[0]) {
            return;
        }

        // Need C up to required N (pdf uses n+1; derivative uses same)
        ensureCForTime(branchLength, /*extraN=*/1);
        final int N = computedN - 1;

        // xh depends on y with dxh/dy = invRewardDiff[h]
        final double dxh_dy = invRewardDiff[h];

        // Poisson / uniformization terms
        final double lt = lambda * branchLength;
        double premult = Math.exp(-lt); // n=0 premult

        // reuse inc buffer (SORTED uv)
        final double[] incSorted = s.inc;

        // Bernstein recurrence helpers
        final double xh = s.xh[0];
        final double oneMinus = 1.0 - xh;
        final double ratio = xh / oneMinus;

        // Track (1-xh)^(n-1) for degree (n-1) Bernstein polynomials:
        // for n=1, (1-xh)^(0)=1
        double w0m = 1.0;

        for (int n = 0; n <= N; ++n) {
            if (n >= 1) {
                // d/dy f_Y = [t * (lambda * invDiff[h]) * premult] * [dxh/dy] * [d/dxh ...]
                // and d/dxh introduces factor n and second differences in k:
                // scale = t * lambda * premult * n * (invDiff[h]^2)
                final double scale = branchLength * lambda * premult * n * (invRewardDiff[h] * dxh_dy);

                Arrays.fill(incSorted, 0.0);

                // degree (n-1) weights
                double w = w0m;

                for (int k = 0; k <= n - 1; ++k) {
                    // Second forward diff of C at (n+1, ·):
                    // C(k+2) - 2 C(k+1) + C(k)
                    final int c0 = cOffset(h, n + 1, k);
                    final int c1 = cOffset(h, n + 1, k + 1);
                    final int c2 = cOffset(h, n + 1, k + 2);

                    for (int uv = 0; uv < dim2; ++uv) {
                        incSorted[uv] += w * (Cflat[c2 + uv] - 2.0 * Cflat[c1 + uv] + Cflat[c0 + uv]);
                    }

                    // advance w -> w_{k+1} for degree (n-1)
                    if (k < n - 1) {
                        w *= ((double) ((n - 1) - k) / (double) (k + 1)) * ratio;
                    }
                }

                // write to ORIGINAL order
                for (int uS = 0; uS < dim; ++uS) {
                    final int outRowBase = outRowBaseBySorted[uS];
                    final int inRowBase = uS * dim;
                    for (int vS = 0; vS < dim; ++vS) {
                        out[outRowBase + outColBySorted[vS]] += scale * incSorted[inRowBase + vS];
                    }
                }
            }

            // premult_{n+1} = premult_n * (lambda*t)/(n+1)
            premult *= lt / (n + 1.0);

            // update w0m when moving n -> n+1:
            // exponent (n) - 1 increases by 1 for n>=1
            if (n >= 1) {
                w0m *= oneMinus;
            }
        }
        if (dim == 2) {
            final double[] fY = computePdf(y, branchLength);
            final double det = fY[0]*fY[3] - fY[1]*fY[2];
            final double d00 = out[0], d01 = out[1], d10 = out[2], d11 = out[3];
            out[0] = ( fY[3]*d00 - fY[1]*d10) / det;   // was fY[2], must be fY[1]
            out[1] = ( fY[3]*d01 - fY[1]*d11) / det;   // was fY[2], must be fY[1]
            out[2] = (-fY[2]*d00 + fY[0]*d10) / det;   // was fY[1], must be fY[2]
            out[3] = (-fY[2]*d01 + fY[0]*d11) / det;   // was fY[1], must be fY[2]
        }
//        if (dim == 2) {
//            // Get the PDF matrix f_Y(y|t) — this is what BEAGLE has as the transition matrix
//            final double[] fY = computePdf(y, branchLength);
//            // 2x2 solve: out = f_Y(y|t)^{-1} @ out  via Cramer's rule
//            final double det = fY[0]*fY[3] - fY[1]*fY[2];
//            final double d00 = out[0], d01 = out[1], d10 = out[2], d11 = out[3];
//            out[0] = ( fY[3]*d00 - fY[2]*d10) / det;
//            out[1] = ( fY[3]*d01 - fY[2]*d11) / det;
//            out[2] = (-fY[1]*d00 + fY[0]*d10) / det;
//            out[3] = (-fY[1]*d01 + fY[0]*d11) / det;
//        }
    }
    /**
     * Compute the "PDF inner increment" in SORTED uv space:
     *
     *   incSorted[uv] = Σ_{k=0}^n B_{n,k}(x_h) * ( C(h,n+1,k+1,uv) - C(h,n+1,k,uv) )
     *
     * using the same Bernstein recurrence as loopCyclePdfAddToOriginalOrderFast_Y,
     * but WITHOUT applying any prefactor and WITHOUT writing to ORIGINAL order.
     */
    private void computePdfInnerIncSorted_Y(
            int h, int n,
            double ratio,
            double w0,            // (1-xh)^n
            double[] incSorted) {

        Arrays.fill(incSorted, 0.0);

        // k = 0
        {
            final int aOff = cOffset(h, n + 1, 1);
            final int bOff = cOffset(h, n + 1, 0);
            for (int uv = 0; uv < dim2; ++uv) {
                incSorted[uv] += w0 * (Cflat[aOff + uv] - Cflat[bOff + uv]);
            }
        }

        // k = 1..n via Bernstein recurrence
        double w = w0;
        for (int k = 0; k < n; ++k) {
            w *= ((double) (n - k) / (double) (k + 1)) * ratio;
            final int kp1 = k + 1;

            final int aOff = cOffset(h, n + 1, kp1 + 1);
            final int bOff = cOffset(h, n + 1, kp1);

            for (int uv = 0; uv < dim2; ++uv) {
                incSorted[uv] += w * (Cflat[aOff + uv] - Cflat[bOff + uv]);
            }
        }
    }
    /**
     * Derivative of the BRANCH-NORMALIZED reward PDF matrix w.r.t. time t (branchLength):
     *
     *   Computes d/dt f_{ij}(y | t) where y = totalReward / t.
     *
     * Matches LaTeX Eq. (forwardGradientWrtBranchLength) with:
     *   - premult = p_n(t) = e^{-λt} (λt)^n / n!
     *   - dp_n/dt = λ (p_{n-1}(t) - p_n(t)), with p_{-1}=0
     *   - Δ_h = r_h - r_{h-1}, and invRewardDiff[h] = 1/Δ_h
     *
     * IMPORTANT:
     *   This assumes h is fixed for the given y (i.e., y not moving across knots).
     *   If y is exactly at a knot, getHfromY(y) chooses the upper interval (your current rule).
     */
    public void computePdfDerivativeWrtTimeInto(double y, double branchLength, double[] out) {
        if (out == null || out.length != dim2) {
            throw new IllegalArgumentException("out must be length dim*dim=" + dim2);
        }
        if (branchLength <= 0.0) {
            throw new IllegalArgumentException("branchLength must be > 0");
        }
        Arrays.fill(out, 0.0);

        ensureNumericsUpToDate();

        // Determine h from y and precompute xh
        final int h = getHfromY(y);
        final Scratch s = tls.get();
        ensureScratchCapacity(s, 1);
        fillXhPrecompY(s, 0, y, h);

        // Ensure C cache (pdf uses n+1)
        ensureCForTime(branchLength, /*extraN=*/1);
        final int N = computedN - 1;

        final double invDiff = invRewardDiff[h];           // 1/Δ_h
        final double lt = lambda * branchLength;

        // premult = p_n(t), start at n=0
        double premult = Math.exp(-lt);   // p_0(t)
        double prevPremult = 0.0;         // p_{-1}(t) := 0

        // Reuse inc buffer (SORTED uv)
        final double[] incSorted = s.inc;

        // For interior Bernstein recurrence
        final boolean isBoundary = (s.isZero[0] || s.isOne[0]);
        final double xh = s.xh[0];
        final double oneMinus = 1.0 - xh;
        final double ratio = (!isBoundary ? (xh / oneMinus) : 0.0);

        // Track w0 = (1-xh)^n for degree-n Bernstein weights (interior only)
        double w0 = 1.0;

        for (int n = 0; n <= N; ++n) {

            final double deltaPremult = prevPremult - premult; // (p_{n-1} - p_n)

            // Term 1 prefactor: (λ/Δ_h) * p_n(t)
            final double tempA = (lambda * invDiff) * premult;

            // Term 2 prefactor: (λ^2 t / Δ_h) * (p_{n-1}(t) - p_n(t))
            final double tempB = (lambda * invDiff) * (lambda * branchLength) * deltaPremult;

            if (isBoundary) {
                // At xh=0 or xh=1, the Bernstein mass is at k=0 or k=n.
                if (s.isZero[0]) {
                    // k=0 contribution: ΔC = C(n+1,1) - C(n+1,0)
                    addDiffBlockToOriginalOrder(out, tempA, cOffset(h, n + 1, 1), cOffset(h, n + 1, 0));
                    addDiffBlockToOriginalOrder(out, tempB, cOffset(h, n + 1, 1), cOffset(h, n + 1, 0));
                } else {
                    // xh=1: k=n contribution: ΔC = C(n+1,n+1) - C(n+1,n)
                    addDiffBlockToOriginalOrder(out, tempA, cOffset(h, n + 1, n + 1), cOffset(h, n + 1, n));
                    addDiffBlockToOriginalOrder(out, tempB, cOffset(h, n + 1, n + 1), cOffset(h, n + 1, n));
                }
            } else {
                // Interior: compute shared inner sum once, then apply both prefactors
                computePdfInnerIncSorted_Y(h, n, ratio, w0, incSorted);

                for (int uS = 0; uS < dim; ++uS) {
                    final int outRowBase = outRowBaseBySorted[uS];
                    final int inRowBase = uS * dim;
                    for (int vS = 0; vS < dim; ++vS) {
                        final int outIdx = outRowBase + outColBySorted[vS];
                        out[outIdx] += (tempA + tempB) * incSorted[inRowBase + vS];
                    }
                }

                // update w0: (1-xh)^n -> (1-xh)^(n+1)
                w0 *= oneMinus;
            }

            // advance Poisson weights: p_{n+1} = p_n * (λt)/(n+1) = p_n * lt/(n+1)
            prevPremult = premult;
            premult *= lt / (n + 1.0);
        }
    }
}