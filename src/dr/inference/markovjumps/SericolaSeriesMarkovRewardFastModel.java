/*
 * SericolaSeriesMarkovRewardFastModel.java
 *
 * Continuous reward density ONLY, in proportion space:
 *   alpha_i in [0,1] (reward-rate proportions)
 *   rho     in (0,1) (reward proportion)
 *
 * Physical mapping (handled outside this class):
 *   a_i = L + w * alpha_i
 *   r   = L + w * rho
 *
 * Assumptions:
 *  - Assumes alpha values are (after sorting) nondecreasing, and strictly increasing
 *    where intervals are used (ties are not allowed where h requires 1/(alpha_h-alpha_{h-1})).
 *  - Endpoints alpha_{i_L}=0 and alpha_{i_U}=1 are allowed.
 *  - rho is expected in [min alpha, max alpha]
 *  - This returns ONLY the continuous component; the atomic mass is handled elsewhere.
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
import dr.math.matrixAlgebra.Vector;

import java.util.Arrays;

/**
 * @author Filippo Monti
 * @author Marc Suchard
 */

public class SericolaSeriesMarkovRewardFastModel extends AbstractModel {

    private static final boolean DEBUG = false;

    // Dependencies
    private final SubstitutionModel underlyingSubstitutionModel;

    // Dimensions
    private final int dim;
    private final int dim2;
    private final int phi; // dim - 1

    // Accuracy / truncation
    private final double epsilon;

    // -----------------------
    // Sorting buffers / maps
    // -----------------------
    private final double[] unsortedQ;     // dim2
    private final double[] sortedQ;       // dim2
    private final double[] sortedAlpha;   // dim

    // perm maps SORTED index -> ORIGINAL index
    private final int[] perm;
    // invPerm maps ORIGINAL index -> SORTED index
    private final int[] invPerm;

    // Mapping helpers for writing results directly to ORIGINAL order
    private final int[] outRowBaseBySorted; // perm[uS] * dim
    private final int[] outColBySorted;     // perm[vS]

    private boolean permDirty = true;
    private boolean qDirty = true;
    private boolean cachesDirty = true;

    // -----------------------
    // Numerics (mutable)
    // -----------------------
    private double lambda;                 // uniformization rate
    private final double[] P;              // dim2, P = I + sortedQ/lambda
    // invAlphaDiff[h] = 1/(sortedAlpha[h]-sortedAlpha[h-1]) for h=1..phi
    private final double[] invAlphaDiff;

    // Exponential for conditional probabilities (optional)
    private final EigenSystem eigenSystem;
    private EigenDecomposition eigenDecomposition;

    // Scratch for computeChnk Pn multiplication
    private final double[] PnScratch;      // dim2
    private final double[] matMulScratch;  // dim2

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

    private int allocatedN = -1;
    private int allocatedN1 = 0;   // allocatedN + 1

    private int computedN = -1;    // -1 means nothing computed / invalid
    private double maxTime = 0.0;  // max time covered by current computed cache
    private final int[] idx;   // length dim, allocated in constructor
    private final double[] out;

    private final boolean conditionalOnZ0;

    // Thread-local scratch to avoid allocations in hot path
    private static final class Scratch {
        int[] H;
        int[] NN;
        double[] lt;
        double[] premult;
        double[] xh;        // in [0,1]
        double[] oneMinus;  // 1-xh
        double[] ratio;     // xh/(1-xh)
        double[] w0;        // (1-xh)^n evolving across n
        boolean[] isZero;   // xh ~ 0
        boolean[] isOne;    // xh ~ 1
        double[] inc;       // length dim2, SORTED uv
    }

    private final ThreadLocal<Scratch> tls = ThreadLocal.withInitial(Scratch::new);

    // Temporary scalar wrappers (allocation-free)
    private final double[] scratchRho1 = new double[1];
    private final double[] scratchT1 = new double[1];
    private final double[][] scratchW1 = new double[1][];

    private final Parameter rewardRatesValues;   // distinct reward values
    private final Parameter rewardRatesValuesInternal;
    private final Parameter rewardRatesMapping;  // state i -> index into rewardRatesValues
    private final double[] alphaVals;

    public SericolaSeriesMarkovRewardFastModel(
            SubstitutionModel underlyingSubstitutionModel,
            Parameter rewardRatesValues,
            Parameter rewardRatesValuesInternal,
            Parameter rewardRatesMapping,
            int dim,
            double epsilon,
            boolean conditionalOnZ0
    ) {
        super("SericolaSeriesMarkovRewardFastModel");

        if (underlyingSubstitutionModel == null) {
            throw new IllegalArgumentException("underlyingSubstitutionModel must be non-null");
        }
        if (rewardRatesValues == null) {
            throw new IllegalArgumentException("rewardRatesValues must be non-null");
        }
        if (rewardRatesMapping == null) {
            throw new IllegalArgumentException("rewardRatesMapping must be non-null");
        }
        if (dim <= 0) {
            throw new IllegalArgumentException("dim must be > 0");
        }
        if (!(epsilon > 0.0 && epsilon < 1.0)) {
            throw new IllegalArgumentException("epsilon must be in (0,1)");
        }

        this.underlyingSubstitutionModel = underlyingSubstitutionModel;
        this.rewardRatesValues = rewardRatesValues;
        this.rewardRatesValuesInternal = rewardRatesValuesInternal;
        this.rewardRatesMapping = rewardRatesMapping;
        this.alphaVals = new double[dim];
        for (int i = 0; i < rewardRatesValues.getDimension(); i++) {
            final double v = rewardRatesValues.getParameterValue(i);
            if (v < 0 || v > 1) {
                throw new IllegalArgumentException("rewardRatesValues must be in [0,1], but value at index " + i + " is " + v);
            }
        }

        this.dim = dim;
        this.dim2 = dim * dim;
        this.phi = dim - 1;
        this.epsilon = epsilon;

        this.unsortedQ = new double[dim2];
        this.sortedQ = new double[dim2];
        this.sortedAlpha = new double[dim];

        this.perm = new int[dim];
        this.invPerm = new int[dim];
        this.outRowBaseBySorted = new int[dim];
        this.outColBySorted = new int[dim];

        this.P = new double[dim2];
        this.invAlphaDiff = new double[dim];

        this.eigenSystem = new DefaultEigenSystem(dim);
        this.PnScratch = new double[dim2];
        this.matMulScratch = new double[dim2];
        this.out = new double[dim2];
        this.conditionalOnZ0 = conditionalOnZ0; // TODO add this to the signature

        addModel(underlyingSubstitutionModel);
        addVariable(rewardRatesValues);
        addVariable(rewardRatesValuesInternal);
        addVariable(rewardRatesMapping);

        this.permDirty = true;
        this.qDirty = true;
        this.cachesDirty = true;
        this.idx = new int[dim];
        for (int i = 0; i < dim; i++) idx[i] = i;   // do this ONCE

    }

    // ============================================================
    // Public API: continuous density f^*(rho, t)
    // ============================================================

    public double[] computePdf(double rho, double time) {
        computePdfInto(rho, time, out);
        return out;
    }

    public void computePdfInto(double rho, double time, double[] out) {
        if (out == null || out.length != dim2) {
            throw new IllegalArgumentException("out must be length dim*dim=" + dim2);
        }
        scratchRho1[0] = rho;
        scratchT1[0] = time;
        scratchW1[0] = out;
        computePdfInto(scratchRho1, scratchT1, false, scratchW1);
    }

    public void computePdfInto(double[] RHO, double time, double[][] W) {
        computePdfInto(RHO, new double[]{time}, false, W);
    }

    public void computePdfInto(double[] RHO, double[] times, double[][] W) {
        computePdfInto(RHO, times, false, W);
    }

    public void computePdfInto(double[] RHO, double[] times, boolean parsimonious, double[][] W) {
        validateInputs(RHO, times, W);
        final int T = RHO.length;
        if (T == 0) return;

        ensureNumericsUpToDate();
        validateAndZeroOutputs(W, T);

        final Scratch s = tls.get();
        ensureScratchCapacity(s, T);

        final double maxT = precomputePdfScratch(RHO, times, parsimonious, s);

        // Ensure C is available up to required N (pdf uses n+1)
        ensureCForTime(maxT, /*extraN=*/1);
        final int N = computedN - 1;

        accumulatePdfOverN(W, s, T, N);

        if (conditionalOnZ0) {
            if (times.length == 1) {
                final double denom = -Math.expm1(-lambda * times[0]); // 1 - exp(-lambda t)
                if (denom <= 0.0) throw new IllegalStateException("Invalid denominator for conditional probability: " + denom);
                final double inv = 1.0 / denom;
                for (int t = 0; t < T; ++t) {
                    final double[] row = W[t];
                    for (int k = 0; k < dim2; ++k) row[k] *= inv;
                }
            } else {
                for (int t = 0; t < T; ++t) {
                    final double denom = -Math.expm1(-lambda * times[t]);
                    if (denom <= 0.0) throw new IllegalStateException("Invalid denominator for conditional probability at index " + t + ": " + denom);
                    final double inv = 1.0 / denom;
                    final double[] row = W[t];
                    for (int k = 0; k < dim2; ++k) row[k] *= inv;
                }
            }
        }
    }

    private void validateInputs(double[] RHO, double[] times, double[][] W) {
        if (RHO == null || times == null || W == null) {
            throw new IllegalArgumentException("RHO/times/W must be non-null");
        }
        final int T = RHO.length;
        final boolean singleTime = (times.length == 1);
        if (!singleTime && times.length != T) {
            throw new IllegalArgumentException("Either times.length==1 or times.length==RHO.length");
        }
        if (W.length != T) {
            throw new IllegalArgumentException("W.length must equal RHO.length");
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

    private double precomputePdfScratch(double[] RHO, double[] times, boolean parsimonious, Scratch s) {
        final int T = RHO.length;
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

                final double rho = RHO[t];
                final int h = getHfromRho(rho);
                s.H[t] = h;

                fillXhPrecomp(s, t, rho, h);
            }
            return maxT;
        }

        for (int t = 0; t < T; ++t) {
            final double time = times[t];
            if (time <= 0.0) throw new IllegalArgumentException("time must be > 0 at index " + t);
            if (time > maxT) maxT = time;

            final double lt = lambda * time;
            s.lt[t] = lt;
            s.premult[t] = Math.exp(-lt);
            s.NN[t] = (parsimonious ? determineNumberOfSteps(time) : Integer.MAX_VALUE);

            final double rho = RHO[t];
            final int h = getHfromRho(rho);
            s.H[t] = h;

            fillXhPrecomp(s, t, rho, h);
        }
        return maxT;
    }

    private void fillXhPrecomp(Scratch s, int t, double rho, int h) {
        final double invDiff = invAlphaDiff[h];
        double xh = (rho - sortedAlpha[h - 1]) * invDiff;

        if (xh <= 0.0) {
            s.xh[t] = 0.0;
            s.isZero[t] = true;
            s.isOne[t] = false;
            return;
        }
        if (xh >= 1.0) {
            s.xh[t] = 1.0;
            s.isZero[t] = false;
            s.isOne[t] = true;
            return;
        }

        s.xh[t] = xh;
        s.isZero[t] = false;
        s.isOne[t] = false;

        final double om = 1.0 - xh;
        s.oneMinus[t] = om;
        s.ratio[t] = xh / om;
        s.w0[t] = 1.0;
    }

    private void accumulatePdfOverN(double[][] W, Scratch s, int T, int N) {
        for (int n = 0; n <= N; ++n) {

            for (int t = 0; t < T; ++t) {
                if (s.NN[t] < n) continue;

                final int h = s.H[t];

                // f^*(rho,t): prefactor is (lambda * time) / (alpha_h - alpha_{h-1}) * p_n(t)
                final double time = s.lt[t] / lambda;
                final double base = (lambda * invAlphaDiff[h]) * s.premult[t] * time;

                if (s.isZero[t]) {
                    addDiffBlockToOriginalOrder(W[t], base,
                            cOffset(h, n + 1, 1),
                            cOffset(h, n + 1, 0));
                } else if (s.isOne[t]) {
                    addDiffBlockToOriginalOrder(W[t], base,
                            cOffset(h, n + 1, n + 1),
                            cOffset(h, n + 1, n));
                } else {
                    loopCyclePdfAddToOriginalOrderFast(
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

    private void loopCyclePdfAddToOriginalOrderFast(
            int h, int n,
            double ratio,
            double w0,                 // (1-xh)^n
            double premult,
            double time,
            double[] WtOriginal,
            double[] incSorted) {

        final double temp = (lambda * invAlphaDiff[h]) * premult * time;

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

    private void addDiffBlockToOriginalOrder(double[] WtOriginal, double temp, int aOff, int bOff) {
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

    // ============================================================
    // C cache growth + computation (same recursions; now uses alpha)
    // ============================================================

    private void ensureCForTime(double time, int extraN) {
        if (time <= 0.0) throw new IllegalArgumentException("time must be > 0");

        final int requiredN = determineNumberOfSteps(time) + extraN;

        ensureCapacityN(requiredN);

        if (time > maxTime) maxTime = time;

        if (computedN < requiredN) {
            computeChnk(requiredN);
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

        computedN = -1;
        maxTime = 0.0;
        invalidateDerivativeCache();
    }

    private void computeChnk(int N) {
        if (Cflat == null || allocatedN < 0) {
            throw new IllegalStateException("C storage not initialized");
        }
        clearCflatUpToN(N);

        Arrays.fill(PnScratch, 0.0);
        for (int u = 0; u < dim; ++u) {
            PnScratch[idx(u, u)] = 1.0;
        }

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

                        final double cScalar = (sortedAlpha[u] - sortedAlpha[h]) / (sortedAlpha[u] - sortedAlpha[h - 1]);
                        final double dScalar = (sortedAlpha[h] - sortedAlpha[h - 1]) / (sortedAlpha[u] - sortedAlpha[h - 1]);

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

            rightMultiply(PnScratch, P, matMulScratch);
            System.arraycopy(matMulScratch, 0, PnScratch, 0, dim2);

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

                        final double cScalar = (sortedAlpha[h - 1] - sortedAlpha[u]) / (sortedAlpha[h] - sortedAlpha[u]);
                        final double dScalar = (sortedAlpha[h] - sortedAlpha[h - 1]) / (sortedAlpha[h] - sortedAlpha[u]);

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

        Arrays.fill(invAlphaDiff, 0.0);
        for (int h = 1; h < dim; h++) {
            final double d = sortedAlpha[h] - sortedAlpha[h - 1];
            if (!(d > 0.0)) {
                throw new IllegalArgumentException("alphaRates must be strictly increasing after sorting; tie at h=" + h);
            }
            invAlphaDiff[h] = 1.0 / d;
        }

        eigenDecomposition = null;

        invalidateCComputedExtent();
        cachesDirty = false;
    }

    private void fillStateRewardRates(double[] stateRewardRates) {
        final double[] values = rewardRatesValues.getParameterValues();

        if (rewardRatesMapping.getDimension() != dim) {
            throw new IllegalArgumentException(
                    "rewardRatesMapping length mismatch: expected " + dim +
                            " got " + rewardRatesMapping.getDimension());
        }

        for (int i = 0; i < dim; i++) {
            final double raw = rewardRatesMapping.getParameterValue(i);
            final int map = (int) Math.round(raw);

            if (Math.abs(raw - map) > 1e-10) {
                throw new IllegalArgumentException(
                        "rewardRatesMapping must contain integer indices; found " + raw + " at state " + i);
            }
            if (map < 0 || map >= values.length) {
                throw new IllegalArgumentException(
                        "rewardRatesMapping[" + i + "]=" + map +
                                " is out of bounds for rewardRatesValues of length " + values.length);
            }

            final double a = values[map];
            if (a < 0.0 || a > 1.0 || Double.isNaN(a)) {
                throw new IllegalArgumentException(
                        "rewardRatesValues must lie in [0,1]; found " + a + " at value index " + map);
            }

            stateRewardRates[i] = a;
        }
    }

    private void ensurePermutationUpToDate() {
        if (!permDirty) return;

//        final double[] alphaVals = alphaRates.getParameterValues();
        fillStateRewardRates(alphaVals);

        if (alphaVals.length != dim) {
            throw new IllegalArgumentException("alphaRates length mismatch: expected " + dim + " got " + alphaVals.length);
        }

        for (int i = 0; i < dim; i++) {
            final double a = alphaVals[i];
            if (a < 0.0 || a > 1.0 || Double.isNaN(a)) {
                throw new IllegalArgumentException("alphaRates must lie in [0,1]; found " + a + " at index " + i);
            }
        }

//        Integer[] idx = new Integer[dim];
//        for (int i = 0; i < dim; i++) idx[i] = i;
//        Arrays.sort(idx, Comparator.comparingDouble(i -> alphaVals[i]));
        // idx already contains 0..dim-1 from the constructor
        sortIndicesByKey(idx, alphaVals);  // custom primitive sort

        for (int s = 0; s < dim; s++) {
            perm[s] = idx[s];
            sortedAlpha[s] = alphaVals[perm[s]];
        }
        for (int o = 0; o < dim; o++) {
            invPerm[o] = -1;
        }
        for (int s = 0; s < dim; s++) {
            invPerm[perm[s]] = s;
        }

        for (int uS = 0; uS < dim; uS++) {
            outRowBaseBySorted[uS] = perm[uS] * dim;
            outColBySorted[uS] = perm[uS];
        }

        permDirty = false;
        qDirty = true;
        cachesDirty = true;

        if (DEBUG) {
//            System.err.println("[DEBUG] alphaRates (original): " + Arrays.toString(alphaRates.getParameterValues()));
            System.err.println("[DEBUG] sortedAlpha: " + Arrays.toString(sortedAlpha));
            System.err.println("[DEBUG] perm (sorted -> original): " + Arrays.toString(perm));
            System.err.println("[DEBUG] invPerm (original -> sorted): " + Arrays.toString(invPerm));
        }
    }
    private static void sortIndicesByKey(int[] idx, double[] key) {
        quickSortByKey(idx, key, 0, idx.length - 1);
    }

    private static void quickSortByKey(int[] a, double[] key, int lo, int hi) {
        while (lo < hi) {
            int i = lo, j = hi;
            final double pivot = key[a[(lo + hi) >>> 1]];

            while (i <= j) {
                while (key[a[i]] < pivot) i++;
                while (key[a[j]] > pivot) j--;
                if (i <= j) {
                    final int tmp = a[i];
                    a[i] = a[j];
                    a[j] = tmp;
                    i++;
                    j--;
                }
            }

            // Recurse on smaller partition first (limits stack depth)
            if (j - lo < hi - i) {
                if (lo < j) quickSortByKey(a, key, lo, j);
                lo = i;
            } else {
                if (i < hi) quickSortByKey(a, key, i, hi);
                hi = j;
            }
        }
    }

    private void ensureSortedQUpToDate() {
        if (!qDirty) return;

        underlyingSubstitutionModel.getInfinitesimalMatrix(unsortedQ);
        if (unsortedQ.length != dim2) {
            throw new IllegalStateException("Internal unsortedQ length mismatch");
        }

        for (int iS = 0; iS < dim; ++iS) {
            final int iO = perm[iS];
            final int rowS = iS * dim;
            final int rowO = iO * dim;
            for (int jS = 0; jS < dim; ++jS) {
                final int jO = perm[jS];
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
        final int end = cOffset(phi, N, N) + dim2;
        Arrays.fill(Cflat, 0, end, 0.0);
    }

    private double determineLambda(double[] QrowMajor) {
        double minDiag = QrowMajor[0];
        for (int i = 1; i < dim; ++i) {
            final double d = QrowMajor[i * dim + i];
            if (d < minDiag) minDiag = d;
        }
        return -minDiag;
    }

    private void invalidateCComputedExtent() {
        computedN = -1;
        maxTime = 0.0;
        invalidateDerivativeCache();
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

    public double getLambda() {
        ensureNumericsUpToDate();
        return lambda;
    }

    public double getUniformizationRate() {
        return getLambda();
    }

    private int getHfromRho(double rho) {
        final double lo = sortedAlpha[0];
        final double hi = sortedAlpha[phi];

        if (rho < lo) {
            System.out.println("rho=" + rho + " < min(alpha)=" + lo + "; snapping to boundary");
            throw new IllegalArgumentException("rho must be >= min(alpha)");
        }
        if (rho > hi) {
            System.out.println("rho=" + rho + " > max(alpha)=" + hi + "; snapping to boundary");
            throw new IllegalArgumentException("rho must be <= max(alpha)");
        }

        int idx = Arrays.binarySearch(sortedAlpha, 0, phi + 1, rho);

        if (idx >= 0) {
            if (idx == 0) return 1;
            if (idx == phi) return phi;
            return idx + 1;
        }

        int insertionPoint = -idx - 1;
        return Math.max(1, Math.min(phi, insertionPoint));
    }

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
    // Optional: conditional probabilities exp(Q t)
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
        if (variable == rewardRatesMapping || variable == rewardRatesValues || variable == rewardRatesValuesInternal) {
            permDirty = true;
            qDirty = true;
            cachesDirty = true;
            fireModelChanged();
        }
    }

    @Override
    protected void storeState() {
        // nothing needed (we recompute on restore)
    }

    @Override
    protected void restoreState() {
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

    @Override
    public String toString() {
        ensureNumericsUpToDate();
        StringBuilder sb = new StringBuilder();
        sb.append("sortedQ: ").append(new Vector(sortedQ)).append("\n");
        sb.append("sortedAlpha: ").append(new Vector(sortedAlpha)).append("\n");
        sb.append("lambda: ").append(lambda).append("\n");
        sb.append("allocatedN: ").append(allocatedN).append("\n");
        sb.append("maxTime: ").append(maxTime).append("\n");
        sb.append("cprob at maxTime: ").append(new Vector(computeConditionalProbabilities(maxTime))).append("\n");
        return sb.toString();
    }

    public void computePdfDerivativeWrtYInto(double RHO, double time, double[] differential) {
        throw new UnsupportedOperationException("Not implemented yet");
    }

    // ============================================================
    // NO-ALLOCATION derivative of reward-density w.r.t. PROPORTION rho
    // ============================================================
    //
    // Computes: d/drho f^*_{ij}(rho | t)  where rho in (0,1)
    //
    // Convention for boundary snapping (xh==0 or 1):
    //   returns 0 derivative (consistent with your "snap-to-boundary" semantics).

    public void computePdfDerivativeWrtRhoInto(double rho, double branchLength, double[] out, boolean shiftback) {
        if (rho<0){
            System.out.println("rho <0");
        }
        if (out == null || out.length != dim2) {
            throw new IllegalArgumentException("out must be length dim*dim=" + dim2);
        }
        if (branchLength <= 0.0) {
            throw new IllegalArgumentException("branchLength must be > 0");
        }
        Arrays.fill(out, 0.0);
        ensureNumericsUpToDate();
        final int h = getHfromRho(rho);
        final Scratch s = tls.get();
        ensureScratchCapacity(s, 1);
        fillXhPrecomp(s, 0, rho, h);
        if (s.isZero[0] || s.isOne[0]) {
//            return;
            System.out.println("Boundary Values for rho touched");
        }

        final boolean xIsZero = s.isZero[0];
        final boolean xIsOne  = s.isOne[0];

        ensureCForTime(branchLength, /*extraN=*/1);
        ensureDerivativeCacheCapacity();
        final int N = computedN - 1;
        final double dxh_drho = invAlphaDiff[h];
        final double lt = lambda * branchLength;
        double premult = Math.exp(-lt); // n=0 premult
        final double[] incSorted = s.inc;
        final double xh = s.xh[0];
        final double oneMinus = 1.0 - xh;
//        final double ratio = xh / oneMinus;
        final double ratio = (!xIsZero && !xIsOne) ? (xh / oneMinus) : 0.0;

        double w0m = 1.0;

        for (int n = 0; n <= N; ++n) {
            if (n >= 1) {
                // f^*(rho,t) has prefactor (lambda * invDiff[h]) * premult * t, with t = branchLength
                // d/drho introduces dxh/drho and a factor n with second differences in k:
                // scale = t * lambda * premult * n * (invDiff[h]^2)
                final double scale = branchLength * lambda * premult * n * (invAlphaDiff[h] * dxh_drho);
                Arrays.fill(incSorted, 0.0);
                if (xIsZero) {
                    // x=0 => only k=0 weight survives, w=1
                    final int k = 0;
                    final int c0 = cOffset(h, n + 1, k);
                    final int c1 = cOffset(h, n + 1, k + 1);
                    final int c2 = cOffset(h, n + 1, k + 2);
                    for (int uv = 0; uv < dim2; ++uv) {
                        incSorted[uv] += (Cflat[c2 + uv] - 2.0 * Cflat[c1 + uv] + Cflat[c0 + uv]);
                    }

                } else if (xIsOne) {
                    // x=1 => only k=n-1 weight survives, w=1
                    final int k = n - 1;
                    final int c0 = cOffset(h, n + 1, k);
                    final int c1 = cOffset(h, n + 1, k + 1);
                    final int c2 = cOffset(h, n + 1, k + 2);
                    for (int uv = 0; uv < dim2; ++uv) {
                        incSorted[uv] += (Cflat[c2 + uv] - 2.0 * Cflat[c1 + uv] + Cflat[c0 + uv]);
                    }
                } else {
                    prepareD2Row(h, n);
                    final int d2Base = d2Offset(h, n, 0);
                    final double[] d2 = D2flat;
                    final double[] inc = incSorted;

                    double w = w0m;
                    for (int k = 0; k <= n - 1; ++k) {
                        final int off = d2Base + k * dim2;
                        final double wk = w;

                        int uv = 0;
                        for (; uv <= dim2 - 4; uv += 4) {
                            inc[uv    ] += wk * d2[off + uv    ];
                            inc[uv + 1] += wk * d2[off + uv + 1];
                            inc[uv + 2] += wk * d2[off + uv + 2];
                            inc[uv + 3] += wk * d2[off + uv + 3];
                        }
                        for (; uv < dim2; ++uv) {
                            inc[uv] += wk * d2[off + uv];
                        }

                        if (k < n - 1) {
                            w *= ((double) (n - 1 - k) / (double) (k + 1)) * ratio;
                        }
                    }
                }

                for (int uS = 0; uS < dim; ++uS) {
                    final int outRowBase = outRowBaseBySorted[uS];
                    final int inRowBase = uS * dim;
                    for (int vS = 0; vS < dim; ++vS) {
                        out[outRowBase + outColBySorted[vS]] += scale * incSorted[inRowBase + vS];
                    }
                }
            }
            premult *= lt / (n + 1.0);
            if (n >= 1) {
                w0m *= oneMinus;
            }
        }
    }

    // Derivative-only lazy cache for second differences in k
    private double[] D2flat;     // stores blocks of length dim2
    private int[] d2Stamp;
    private int d2Epoch = 1;
    private int d2AllocatedN = -1;

    private int d2BlockIndex(int h, int n, int k) {
        return ((h * allocatedN1 + n) * allocatedN1 + k);
    }

    private int d2Offset(int h, int n, int k) {
        return d2BlockIndex(h, n, k) * dim2;
    }
//    private void invalidateDerivativeCache() {
//        D2flat = null;
//        d2Ready = null;
//        d2AllocatedN = -1;
//    }
private void invalidateDerivativeCache() {
    if (d2Stamp == null) {
        return;
    }
    d2Epoch++;
    if (d2Epoch == Integer.MAX_VALUE) {
        Arrays.fill(d2Stamp, 0);
        d2Epoch = 1;
    }
}
//    private void ensureDerivativeCacheCapacity() {
//        if (allocatedN < 0) {
//            throw new IllegalStateException("C cache must be allocated first");
//        }
//
//        if (D2flat != null && d2AllocatedN == allocatedN) {
//            return;
//        }
//
//        final long nBlocks = (long) (phi + 1) * (long) allocatedN1 * (long) allocatedN1;
//        final long totalDoubles = nBlocks * (long) dim2;
//
//        if (totalDoubles > Integer.MAX_VALUE) {
//            throw new IllegalStateException("D2 array too large: " + totalDoubles + " doubles");
//        }
//
//        D2flat = new double[(int) totalDoubles];
//        d2Ready = new byte[(int) nBlocks];
//        d2AllocatedN = allocatedN;
//    }

    private void ensureDerivativeCacheCapacity() {

        if (allocatedN < 0) {
            throw new IllegalStateException("C cache must be allocated first");
        }

        if (D2flat != null && d2AllocatedN == allocatedN) {
            return;
        }

        final long nBlocks = (long) (phi + 1) * (long) allocatedN1 * (long) allocatedN1;
        final long totalDoubles = nBlocks * (long) dim2;

        System.out.println("ALLOC D2: allocatedN=" + allocatedN
                + " allocatedN1=" + allocatedN1
                + " nBlocks=" + nBlocks
                + " totalDoubles=" + totalDoubles
                + " MB=" + (8.0 * totalDoubles) / (1024.0 * 1024.0));

        if (totalDoubles > Integer.MAX_VALUE) {
            throw new IllegalStateException("D2 array too large: " + totalDoubles + " doubles");
        }
        if (nBlocks > Integer.MAX_VALUE) {
            throw new IllegalStateException("d2Stamp array too large: " + nBlocks + " blocks");
        }

        D2flat = new double[(int) totalDoubles];
        d2Stamp = new int[(int) nBlocks];
        d2AllocatedN = allocatedN;
        d2Epoch = 1;
    }


    private void prepareD2Row(int h, int n) {
        final int n1 = n + 1;
        final int baseBlock = d2BlockIndex(h, n, 0);
        final int baseD2Off = baseBlock * dim2;

        for (int k = 0; k <= n - 1; ++k) {
            final int block = baseBlock + k;
            if (d2Stamp[block] != d2Epoch) {
                final int d2Off = baseD2Off + k * dim2;
                final int c0 = cOffset(h, n1, k);
                final int c1 = cOffset(h, n1, k + 1);
                final int c2 = cOffset(h, n1, k + 2);

                final double[] C = Cflat;
                final double[] D2 = D2flat;

                for (int uv = 0; uv < dim2; ++uv) {
                    D2[d2Off + uv] =
                            C[c2 + uv]
                                    - 2.0 * C[c1 + uv]
                                    + C[c0 + uv];
                }
//                int uv = 0;
//                for (; uv <= dim2 - 4; uv += 4) {
//                    D2[d2Off + uv]     = C[c2 + uv]     - 2.0*C[c1 + uv]     + C[c0 + uv];
//                    D2[d2Off + uv + 1] = C[c2 + uv + 1] - 2.0*C[c1 + uv + 1] + C[c0 + uv + 1];
//                    D2[d2Off + uv + 2] = C[c2 + uv + 2] - 2.0*C[c1 + uv + 2] + C[c0 + uv + 2];
//                    D2[d2Off + uv + 3] = C[c2 + uv + 3] - 2.0*C[c1 + uv + 3] + C[c0 + uv + 3];
//                }
//                for (; uv < dim2; ++uv) {
//                    D2[d2Off + uv] = C[c2 + uv] - 2.0*C[c1 + uv] + C[c0 + uv];
//                }
                d2Stamp[block] = d2Epoch;
            }
        }
    }


    private int ensureD2Block(int h, int n, int k) {
//        ensureDerivativeCacheCapacity();

        final int block = d2BlockIndex(h, n, k);
        final int d2Off = block * dim2;

        if (d2Stamp[block] != d2Epoch) {
            final int c0 = cOffset(h, n + 1, k);
            final int c1 = cOffset(h, n + 1, k + 1);
            final int c2 = cOffset(h, n + 1, k + 2);

            for (int uv = 0; uv < dim2; ++uv) {
                D2flat[d2Off + uv] =
                        Cflat[c2 + uv]
                                - 2.0 * Cflat[c1 + uv]
                                + Cflat[c0 + uv];
            }
            d2Stamp[block] = d2Epoch;
        }

        return d2Off;
    }



    /**
     * Compute the "PDF inner increment" in SORTED uv space:
     *
     *   incSorted[uv] = Σ_{k=0}^n B_{n,k}(x_h) * ( C(h,n+1,k+1,uv) - C(h,n+1,k,uv) )
     *
     * using the same Bernstein recurrence as loopCyclePdfAddToOriginalOrderFast,
     * but WITHOUT applying any prefactor and WITHOUT writing to ORIGINAL order.
     */
    private void computePdfInnerIncSorted_Rho(
            int h, int n,
            double ratio,
            double w0,            // (1-xh)^n
            double[] incSorted) {
        Arrays.fill(incSorted, 0.0);
        {
            final int aOff = cOffset(h, n + 1, 1);
            final int bOff = cOffset(h, n + 1, 0);
            for (int uv = 0; uv < dim2; ++uv) {
                incSorted[uv] += w0 * (Cflat[aOff + uv] - Cflat[bOff + uv]);
            }
        }
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
}
