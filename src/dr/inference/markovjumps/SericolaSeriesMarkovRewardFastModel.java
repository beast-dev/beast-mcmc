/*
 * SericolaSeriesMarkovRewardFastModel.java
 *
 * A BEAST Model wrapper around the fast Sericola-series Markov reward PDF/CDF engine.
 *
 * Key design:
 *  - Listens to:
 *      (i) underlyingSubstitutionModel  -> provides UNSORTED infinitesimal matrix Q
 *      (ii) rewardRates Parameter       -> provides UNSORTED reward rates r
 *  - Internally computes permutation perm that sorts rewards increasingly, then builds:
 *      sortedR, sortedQ, lambda, P, invRewardDiff, (optional eigenDecomposition)
 *  - Cflat cache is computed in SORTED index space.
 *  - Public outputs are written in ORIGINAL (unsorted) state order WITHOUT allocating
 *    per-branch temporary matrices: we map sorted (uS,vS) -> original (uO,vO) on the fly.
 *
 * Store/restore:
 *  - We cannot snapshot Cflat (potentially huge).
 *  - On restore we "invalidate computed extent" but keep allocations to avoid reallocation.
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
    private int capacityN = -1; // maximum n currently computed (>=0), -1 means none computed
    private int N1 = 0;         // = capacityN + 1

    // Tracks max time seen in current computed cache extent
    private double maxTime = 0.0;

    // Thread-safe scratch buffers to avoid allocations in computePdf hot path
    private static final class Scratch {
        int[] H;
        int[] NN;
        double[] lt;
        double[] premult;

        // NEW: precomputed per-t
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
    private final double[] scratchX1 = new double[1];
    private final double[] scratchT1 = new double[1];
    private final double[][] scratchW1 = new double[1][];

    // ---------------------------------------
    // MCMC store/restore bookkeeping
    // ---------------------------------------
    private boolean storedPermDirty, storedQDirty, storedCachesDirty;
    private int storedCapacityN, storedN1;
    private double storedMaxTime;
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
    // Public API: PDF
    // ============================================================

    @Override
    public double[] computePdf(double reward, double branchLength) {
        final double[] out = new double[dim2];
        computePdfInto(reward, branchLength, out);
        return out;
    }

    @Override
    public double computePdf(double reward, double branchLength, int i, int j) {
        if (i < 0 || i >= dim || j < 0 || j >= dim) {
            throw new IndexOutOfBoundsException("i,j out of bounds");
        }
        final double[] out = new double[dim2];
        computePdfInto(reward, branchLength, out);
        return out[i * dim + j];
    }

    public void computePdfInto(double reward, double branchLength, double[] out) {
        if (out == null || out.length != dim2) {
            throw new IllegalArgumentException("out must be length dim*dim=" + dim2);
        }
        scratchX1[0] = reward;
        scratchT1[0] = branchLength;
        scratchW1[0] = out;
        computePdfInto(scratchX1, scratchT1, false, scratchW1);
    }

    public void computePdfInto(double[] X, double time, double[][] W) {
        computePdfInto(X, new double[]{time}, false, W);
    }

    public void computePdfInto(double[] X, double[] times, double[][] W) {
        computePdfInto(X, times, false, W);
    }

    public void computePdfInto(double[] X, double[] times, boolean parsimonious, double[][] W) {
        if (X == null || times == null || W == null) {
            throw new IllegalArgumentException("X/times/W must be non-null");
        }
        final int T = X.length;
        if (T == 0) return;

        final boolean singleTime = (times.length == 1);
        if (!singleTime && times.length != T) {
            throw new IllegalArgumentException("Either times.length==1 or times.length==X.length");
        }
        if (W.length != T) {
            throw new IllegalArgumentException("W.length must equal X.length");
        }

        // Ensure internal numerics and ordering are up-to-date
        ensureNumericsUpToDate();

        // Validate output rows and zero them (overwrite semantics)
        for (int t = 0; t < T; ++t) {
            final double[] row = W[t];
            if (row == null || row.length != dim2) {
                throw new IllegalArgumentException("W[" + t + "] must be non-null and length dim*dim=" + dim2);
            }
            Arrays.fill(row, 0.0);
        }

        // Thread-local scratch
        final Scratch s = tls.get();
        ensureScratchCapacity(s, T);

        // ---------------------------
        // Precompute per-t invariants
        // ---------------------------
        double maxT = 0.0;

        if (singleTime) {
            final double time = times[0];
            if (time <= 0.0) throw new IllegalArgumentException("time must be > 0");
            maxT = time;

            final double lt = lambda * time;
            final double premult0 = Math.exp(-lt);

            // If you ever want parsimonious in singleTime mode, compute NN0 here.
            final int NN0 = (parsimonious ? determineNumberOfSteps(time) : Integer.MAX_VALUE);

            for (int t = 0; t < T; ++t) {
                final int h = getHfromXSortedRate(X[t], time);
                s.H[t] = h;

                s.lt[t] = lt;
                s.premult[t] = premult0;
                s.NN[t] = NN0;

                // xh in [0,1] + boundary flags + ratio + w0 init
                final double invDiff = invRewardDiff[h];
                double xh = (X[t] - sortedR[h - 1] * time) * invDiff / time;
                if (xh <= 0.0) {
                    xh = 0.0;
                    s.isZero[t] = true;
                    s.isOne[t] = false;
                    s.oneMinus[t] = 1.0;   // unused
                    s.ratio[t] = 0.0;      // unused
                    s.w0[t] = 0.0;         // unused
                } else if (xh >= 1.0) {
                    xh = 1.0;
                    s.isZero[t] = false;
                    s.isOne[t] = true;
                    s.oneMinus[t] = 0.0;   // unused
                    s.ratio[t] = 0.0;      // unused
                    s.w0[t] = 0.0;         // unused
                } else {
                    s.isZero[t] = false;
                    s.isOne[t] = false;
                    final double om = 1.0 - xh;
                    s.oneMinus[t] = om;
                    s.ratio[t] = xh / om;
                    s.w0[t] = 1.0;         // (1-xh)^0
                }

                s.xh[t] = xh; // (optional; mostly for debugging)
            }
        } else {
            for (int t = 0; t < T; ++t) {
                final double time = times[t];
                if (time <= 0.0) throw new IllegalArgumentException("time must be > 0 at index " + t);
                if (time > maxT) maxT = time;

                final int h = getHfromXSortedRate(X[t], time);
                s.H[t] = h;

                final double lt = lambda * time;
                s.lt[t] = lt;
                s.premult[t] = Math.exp(-lt);

                s.NN[t] = (parsimonious ? determineNumberOfSteps(time) : Integer.MAX_VALUE);

                final double invDiff = invRewardDiff[h];
                double xh = (X[t] - sortedR[h - 1] * time) * invDiff / time;
                if (xh <= 0.0) {
                    xh = 0.0;
                    s.isZero[t] = true;
                    s.isOne[t] = false;
                    s.oneMinus[t] = 1.0;
                    s.ratio[t] = 0.0;
                    s.w0[t] = 0.0;
                } else if (xh >= 1.0) {
                    xh = 1.0;
                    s.isZero[t] = false;
                    s.isOne[t] = true;
                    s.oneMinus[t] = 0.0;
                    s.ratio[t] = 0.0;
                    s.w0[t] = 0.0;
                } else {
                    s.isZero[t] = false;
                    s.isOne[t] = false;
                    final double om = 1.0 - xh;
                    s.oneMinus[t] = om;
                    s.ratio[t] = xh / om;
                    s.w0[t] = 1.0; // (1-xh)^0
                }

                s.xh[t] = xh; // (optional)
            }
        }

        // Ensure C is available up to required N (pdf uses n+1)
        ensureCForTime(maxT, /*extraN=*/1);
        final int N = capacityN - 1;

        // ---------------------------
        // Main accumulation over n
        // ---------------------------
        for (int n = 0; n <= N; ++n) {
            for (int t = 0; t < T; ++t) {
                if (s.NN[t] < n) continue;

                final int h = s.H[t];
                final double base = (lambda * invRewardDiff[h]) * s.premult[t];

                if (s.isZero[t]) {
                    final int aOff = cOffset(h, n + 1, 1);
                    final int bOff = cOffset(h, n + 1, 0);
                    addDiffBlockToOriginalOrder(W[t], base, aOff, bOff);
                } else if (s.isOne[t]) {
                    final int aOff = cOffset(h, n + 1, n + 1);
                    final int bOff = cOffset(h, n + 1, n);
                    addDiffBlockToOriginalOrder(W[t], base, aOff, bOff);
                } else {
                    // w0(t) is (1-xh)^n maintained across n; ratio(t) fixed
                    loopCyclePdfAddToOriginalOrderFast(
                            h, n,
                            s.ratio[t],
                            s.w0[t],
                            s.premult[t],
                            W[t],
                            s.inc
                    );
                }
            }

            // premult_{n+1}(t) = premult_n(t) * (lambda*t)/(n+1)
            final double inv = 1.0 / (n + 1.0);
            for (int t = 0; t < T; ++t) {
                s.premult[t] *= s.lt[t] * inv;
            }

            // update w0: (1-xh)^n -> (1-xh)^(n+1)
            for (int t = 0; t < T; ++t) {
                if (!s.isZero[t] && !s.isOne[t]) {
                    s.w0[t] *= s.oneMinus[t];
                }
            }
        }
    }


    private void loopCyclePdfAddToOriginalOrderFast(
            int h, int n,
            double ratio,
            double w0,                 // (1-xh)^n
            double premult,
            double[] WtOriginal,
            double[] incSorted) {

        final double temp = (lambda * invRewardDiff[h]) * premult;

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
    // Public API: CDF (kept; not heavily optimized)
    // ============================================================

    @Override
    public double computeCdf(double x, double time, int i, int j) {
        return computeCdf(new double[]{x}, time)[0][i * dim + j];
    }

    public double[] computeCdf(double x, double time) {
        return computeCdf(new double[]{x}, time)[0];
    }

    public double[][] computeCdf(double[] X, double time) {
        if (time <= 0.0) throw new IllegalArgumentException("time must be > 0");

        ensureNumericsUpToDate();
        ensureCForTime(time, /*extraN=*/0);
        final int N = capacityN;

        final int T = X.length;
        final double[][] W = new double[T][dim2];
        final int[] H = new int[T];
        for (int t = 0; t < T; t++) H[t] = getHfromXSortedRate(X[t], time);

        for (int n = 0; n <= N; ++n) {
            accumulateCdfAddToOriginalOrder(W, X, H, n, time);
        }
        return W;
    }

    // ============================================================
    // Inner kernels
    // ============================================================

    /**
     * Hot inner kernel for PDF: add contribution for given (x,time,h,n) directly into Wt (ORIGINAL order).
     *
     * The cache Cflat is indexed by SORTED uv, but we write into ORIGINAL uv via precomputed maps:
     *   sorted (uS,vS) -> original index (perm[uS], perm[vS]).
     */
    private void loopCyclePdfAddToOriginalOrder(double x, double time, int h, int n,
                                                double premult, double[] WtOriginal, double[] incSorted) {

        final double invDiff = invRewardDiff[h];
        final double factor = lambda * invDiff;

        double xh = (x - sortedR[h - 1] * time) * invDiff / time;
        if (xh <= 0.0) xh = 0.0;
        else if (xh >= 1.0) xh = 1.0;

        final double temp = factor * premult;

        // Boundary cases: avoid Bernstein loop
        if (xh == 0.0) {
            final int aOff = cOffset(h, n + 1, 1);
            final int bOff = cOffset(h, n + 1, 0);
            addDiffBlockToOriginalOrder(WtOriginal, temp, aOff, bOff);
            return;
        }

        if (xh == 1.0) {
            final int aOff = cOffset(h, n + 1, n + 1);
            final int bOff = cOffset(h, n + 1, n);
            addDiffBlockToOriginalOrder(WtOriginal, temp, aOff, bOff);
            return;
        }

        final double oneMinus = 1.0 - xh;
        final double ratio = xh / oneMinus;

        // w0 = (1-x)^n
        double w;
        switch (n) {
            case 0: w = 1.0; break;
            case 1: w = oneMinus; break;
            case 2: w = oneMinus * oneMinus; break;
            case 3: w = oneMinus * oneMinus * oneMinus; break;
            default: w = Math.exp(n * Math.log1p(-xh)); break;
//            default: w = Math.pow(oneMinus, n);
        }

        Arrays.fill(incSorted, 0.0);

        // k=0
        {
            final int aOff = cOffset(h, n + 1, 1);
            final int bOff = cOffset(h, n + 1, 0);
            final int baseA = aOff;
            final int baseB = bOff;
            for (int uv = 0; uv < dim2; ++uv) {
                incSorted[uv] += w * (Cflat[baseA + uv] - Cflat[baseB + uv]);
            }
        }

        // k=1..n via Bernstein recurrence
        for (int k = 0; k < n; ++k) {
            w *= ((double) (n - k) / (double) (k + 1)) * ratio;
            final int kp1 = k + 1;

            final int aOff = cOffset(h, n + 1, kp1 + 1);
            final int bOff = cOffset(h, n + 1, kp1);

            for (int uv = 0; uv < dim2; ++uv) {
                incSorted[uv] += w * (Cflat[aOff + uv] - Cflat[bOff + uv]);
            }
        }

        // Apply temp and write into ORIGINAL order
        // incSorted is in sorted uv: uv = uS*dim + vS
        for (int uS = 0; uS < dim; ++uS) {
            final int outRowBase = outRowBaseBySorted[uS];
            final int inRowBase = uS * dim;
            for (int vS = 0; vS < dim; ++vS) {
                final int outIdx = outRowBase + outColBySorted[vS];
                WtOriginal[outIdx] += temp * incSorted[inRowBase + vS];
            }
        }
    }

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

    // CDF accumulation (structurally similar to your previous version; writes in ORIGINAL order)
    private void accumulateCdfAddToOriginalOrder(double[][] W, double[] X, int[] H, int n, double time) {

        final double lt = lambda * time;
        final double premult = Math.exp(
                -lt + n * (Math.log(lambda) + Math.log(time)) - GammaFunction.lnGamma(n + 1.0)
        );

        for (int t = 0; t < X.length; ++t) {
            final double x = X[t];
            final int h = H[t];

            double xh = (x - sortedR[h - 1] * time) / ((sortedR[h] - sortedR[h - 1]) * time);

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

        if (time > maxTime) maxTime = time;

        final int requiredN = determineNumberOfSteps(time) + extraN;

        if (capacityN < requiredN) {
            ensureCapacityN(requiredN);
            computeChnk(); // recompute all up to capacityN
        }
    }

    private void ensureCapacityN(int requiredN) {
        int newCap = (capacityN < 0) ? requiredN : capacityN;
        while (newCap < requiredN) newCap = Math.max(requiredN, (int)(newCap * 1.5) + 1);

        final int newN1 = newCap + 1;
        final long totalDoubles = (long)(phi + 1) * newN1 * newN1 * dim2;
        if (totalDoubles > Integer.MAX_VALUE) throw new IllegalStateException("C array too large: " + totalDoubles + " doubles"); //TODO THIS SHOULD NOT HAPPEN, block before

        // REUSE if already allocated big enough
        if (Cflat != null && capacityN >= newCap && N1 == newN1) {
            return;
        }

        // If Cflat exists but N1 differs (capacity differs), only allocate when growing
        if (Cflat != null && capacityN >= newCap) {
            // this case shouldn't happen if you keep N1 consistent with capacityN
            return;
        }

        Cflat = new double[(int) totalDoubles];
        capacityN = newCap;
        N1 = newN1;
    }


//    private void ensureCapacityN(int requiredN) {
//        int newCap = (capacityN < 0) ? requiredN : capacityN;
//        while (newCap < requiredN) {
//            newCap = Math.max(requiredN, (int) (newCap * 1.5) + 1);
//        }
//        if (capacityN < 0) newCap = requiredN;
//
//        final int newN1 = newCap + 1;
//
//        final long blocks = (long) (phi + 1) * (long) newN1 * (long) newN1;
//        final long totalDoubles = blocks * (long) dim2;
//
//        if (totalDoubles > Integer.MAX_VALUE) {
//            throw new IllegalStateException("C array too large: " + totalDoubles + " doubles");
//        }
//
//        // Allocate fresh; do NOT copy old (we recompute everything)
//        this.Cflat = new double[(int) totalDoubles];
//        this.capacityN = newCap;
//        this.N1 = newN1;
//
//        if (DEBUG) {
//            System.err.println("Allocated Cflat for capacityN=" + capacityN + " (N1=" + N1 + "), doubles=" + totalDoubles);
//        }
//    }

    /**
     * Compute all C(h,n,k,uv) up to capacityN using Sericola recursions,
     * assuming P, lambda, sortedR are current and consistent.
     *
     * IMPORTANT: Works in SORTED uv space.
     */
    private void computeChnk() {
        if (Cflat == null || capacityN < 0) {
            throw new IllegalStateException("C storage not initialized");
        }

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

        final int N = capacityN;

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

        // Update lambda and P from sortedQ
        this.lambda = determineLambda(sortedQ);
        if (!(lambda > 0.0) || Double.isNaN(lambda) || Double.isInfinite(lambda)) {
            throw new IllegalStateException("Invalid uniformization rate lambda=" + lambda);
        }
        fillP(sortedQ, lambda, P);

        // Update invRewardDiff from sortedR (requires strictly increasing rewards)
        Arrays.fill(invRewardDiff, 0.0);
        for (int h = 1; h < dim; h++) {
            final double d = sortedR[h] - sortedR[h - 1];
            if (!(d > 0.0)) {
                throw new IllegalArgumentException("rewardRates must be strictly increasing after sorting; tie at h=" + h);
            }
            invRewardDiff[h] = 1.0 / d;
        }

        // Invalidate eigen cache
        eigenDecomposition = null;

        // IMPORTANT: invalidate computed extent but keep any allocated Cflat //    todo check the following two lines
//        capacityN = -1;
//        N1 = 0;
        maxTime = 0.0;

        cachesDirty = false;

        if (DEBUG) {
            System.err.println("Sericola numerics refreshed: lambda=" + lambda);
        }
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
    }

    private void ensureSortedQUpToDate() {
        if (!qDirty) return;

        underlyingSubstitutionModel.getInfinitesimalMatrix(unsortedQ);
        if (unsortedQ.length != dim2) {
            throw new IllegalStateException("Internal unsortedQ length mismatch");
        }

        // Reorder Q into sortedQ using perm; matches your previous column-major convention
        for (int j = 0; j < dim; j++) {
            int pj = perm[j];
            for (int i = 0; i < dim; i++) {
                sortedQ[j * dim + i] = unsortedQ[pj * dim + perm[i]];
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
        return (((h * N1 + n) * N1 + k) * dim2);
    }

    private double determineLambda(double[] QrowMajor) {
        double minDiag = QrowMajor[0];          // Q[0,0]
        for (int i = 1; i < dim; ++i) {
            final double d = QrowMajor[i * dim + i];
            if (d < minDiag) minDiag = d;
        }
        return -minDiag;
    }

    private double determineLambda(double[] Qsorted, int dimLocal) {
        double minDiag = Qsorted[0];
        for (int i = 1; i < dimLocal; ++i) {
            final double d = Qsorted[i * dimLocal + i];
            if (d < minDiag) minDiag = d;
        }
        return -minDiag;
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

    /**
     * Determine h based on rate = x/time in SORTED reward space.
     * Returns an h in [1..phi].
     */
    private int getHfromXSortedRate(double x, double time) {
        if (time <= 0.0) throw new IllegalArgumentException("time must be > 0");

        final double lo = sortedR[0] * time;
        final double hi = sortedR[phi] * time;

        if (x < lo) throw new IllegalArgumentException("x must be >= minReward*time");
        if (x > hi) throw new IllegalArgumentException("x must be <= maxReward*time");

        final double rate = x / time;
        int h = Arrays.binarySearch(sortedR, 0, phi + 1, rate);
        if (h < 0) h = -h - 1;

        if (h <= 0) return 1;
        if (h > phi) return phi;
        return h;
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

        storedCapacityN = capacityN;
        storedN1 = N1;
        storedMaxTime = maxTime;

        hasStoredState = true;
    }

    @Override
    protected void restoreState() {
        if (!hasStoredState) return;

        // After BEAST restore, upstream values may have reverted WITHOUT firing events.
        // Therefore, all derived state must be treated as dirty.
        permDirty   = true;
        qDirty      = true;
        cachesDirty = true;

        // Invalidate computed extent, but keep allocations to avoid reallocation.
        // IMPORTANT: do NOT overwrite/assume capacityN is meaningful; it refers to computed extent. // todo check next 3 lines
//        capacityN = -1;
//        N1 = (Cflat == null ? 0 : N1); // keep N1 as-is if you want, but it is irrelevant when capacityN=-1
//        maxTime = 0.0;

        // Eigen cache must be rebuilt too
        eigenDecomposition = null;
    }

//    @Override
//    protected void restoreState() {
//        if (!hasStoredState) return;
//
//        // Restore dirtiness flags, but NEVER trust Cflat contents after restore:
//        // we cannot snapshot it; safest is invalidate computed extent.
//        permDirty = storedPermDirty;
//        qDirty = storedQDirty;
//
//        // Force recompute-on-demand, but keep allocations
//        cachesDirty = true;
//
//        // Invalidate computed extent; keep Cflat allocated as-is
//        capacityN = -1;
//        N1 = 0;
//        maxTime = 0.0;
//    }

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
        sb.append("capacityN: ").append(capacityN).append("\n");
        sb.append("maxTime: ").append(maxTime).append("\n");
        sb.append("cprob at maxTime: ").append(new Vector(computeConditionalProbabilities(maxTime))).append("\n");
        return sb.toString();
    }
}
