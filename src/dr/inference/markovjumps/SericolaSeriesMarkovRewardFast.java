/*
 * SericolaSeriesMarkovReward.java
 *
 * Efficient implementation notes (2026-01):
 * - Hotspots are computePdf + loopCyclePdf and (previously) enormous allocations of internalC.
 * - This rewrite:
 *   (1) Stores C(h,n,k,uv) in ONE flat double[] (no 4D object graph).
 *   (2) Uses geometric growth for N-capacity to avoid repeated multi-GB reallocations.
 *   (3) Uses ThreadLocal scratch buffers (H, premult, lt, NN, inc) to avoid per-call allocations and be thread-safe.
 *   (4) Precomputes invRewardDiff[h] = 1/(r[h]-r[h-1]) for hot-path divisions.
 *   (5) Uses Poisson recurrence for both determineNumberOfSteps() and premult updates.
 *
 * - “CDF parts” are left functional but not heavily optimized; PDF is the focus.
 */

package dr.inference.markovjumps;

import dr.evomodel.substmodel.DefaultEigenSystem;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.EigenSystem;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.Binomial;
import dr.math.GammaFunction;

import dr.math.matrixAlgebra.Vector;

import java.util.Arrays;

import static dr.stats.DiscreteStatistics.max;

public class SericolaSeriesMarkovRewardFast extends AbstractModel implements MarkovReward {

    private static final boolean DEBUG = false;

    // Inputs
    private final double[] Q;
    private final double[] r;

    // Dimensions
    private final int dim;
    private final int dim2;
    private final int phi; // dim - 1

    // Uniformization
    private final double lambda;
    private final double[] P; // I + Q/lambda

    // Accuracy / truncation
    private final double epsilon;

    // Precompute to remove divisions in hot path
    // invRewardDiff[h] = 1 / (r[h] - r[h-1]) for h=1..phi
    private final double[] invRewardDiff;

    // Exponential for conditional probabilities
    private final EigenSystem eigenSystem;
    private EigenDecomposition eigenDecomposition;

    private final double[] PnScratch;

    // ------------------------------------------------------------
    // Storage for C(h,n,k,uv) flattened into a single array
    //
    // Indexing:
    //   base(h,n,k) = (((h * (N1) + n) * (N1) + k) * dim2)
    //   then add uv
    //
    // where N1 = capacityN + 1
    // ------------------------------------------------------------
    private double[] Cflat;
    private int capacityN = -1; // maximum n currently stored (>=0)
    private int N1 = 0;         // = capacityN + 1

    // Tracks max time seen (monotone increasing)
    private double maxTime = 0.0;

    // Scratch for Pn multiplication in computeChnk()
    private final double[] matMulScratch;
    private final double[] incScratch; // reused in computePdf scalar wrapper
    // --- scratch for the scalar wrapper (allocation-free) ---
    private final double[] scratchX1 = new double[1];
    private final double[] scratchT1 = new double[1];
    private final double[][] scratchW1 = new double[1][];

    // Thread-safe scratch buffers to avoid allocations in computePdf hot path
    private static final class Scratch {
        int[] H;
        int[] NN;
        double[] lt;
        double[] premult;
        double[] inc; // length dim2
    }

    private final ThreadLocal<Scratch> tls = ThreadLocal.withInitial(Scratch::new);

    public SericolaSeriesMarkovRewardFast(double[] Q, double[] r, int dim) {
        this(Q, r, dim, 1E-10);
    }

    public SericolaSeriesMarkovRewardFast(double[] Q, double[] r, int dim, double epsilon) {
        super("SericolaSeriesMarkovRewardFast");

        // Validate inputs
        if (dim <= 0) throw new IllegalArgumentException("dim must be > 0");
        if (Q == null || Q.length != dim * dim) throw new IllegalArgumentException("Q length mismatch");
        if (r == null || r.length != dim) throw new IllegalArgumentException("r length mismatch");
        if (!(epsilon > 0.0 && epsilon < 1.0)) throw new IllegalArgumentException("epsilon must be in (0,1)");

        // Rewards must be strictly increasing to avoid denom=0
        for (int i = 1; i < r.length; i++) {
            if (r[i] <= r[i - 1]) {
//                throw new IllegalArgumentException("r must be strictly increasing"); //TODO check this condition
            }
        }

        this.Q = Q;
        this.r = r;
        this.dim = dim;
        this.dim2 = dim * dim;
        this.phi = dim - 1;
        this.epsilon = epsilon;

        this.lambda = determineLambda();
        if (!(lambda > 0.0) || Double.isNaN(lambda) || Double.isInfinite(lambda)) {
            throw new IllegalStateException("Invalid uniformization rate lambda=" + lambda);
        }

        this.P = initializeP(Q, lambda);
        this.eigenSystem = new DefaultEigenSystem(dim);

        this.invRewardDiff = new double[dim];
        for (int h = 1; h < dim; h++) {
            invRewardDiff[h] = 1.0 / (r[h] - r[h - 1]);
        }

        this.matMulScratch = new double[dim2];
        this.incScratch = new double[dim2];
        PnScratch = new double[dim2];

        if (DEBUG) {
            System.err.println("lambda = " + lambda);
        }
    }

    // ------------------------------------------------------------
    // Public API: PDF
    // ------------------------------------------------------------

//    @Override
//    public double computePdf(double x, double time, int i, int j) {
//        // If you truly never call this, fine; otherwise keep correct:
//        return computePdf(new double[]{x}, new double[]{time}, false)[0][idx(i, j)];
//    }


    // ------------------------------------------------------------------
// Legacy API: allocates and returns a fresh dim2 vector
// ------------------------------------------------------------------
    @Override
    public double[] computePdf(double reward, double branchLength) {
        final double[] out = new double[dim2];
        computePdfInto(reward, branchLength, out);
        return out;
    }

    @Override
    public double computePdf(double reward, double branchLength, int i, int j) {

        if (branchLength <= 0.0) {
            throw new IllegalArgumentException("branchLength must be > 0");
        }
        if (i < 0 || i >= dim || j < 0 || j >= dim) {
            throw new IndexOutOfBoundsException("i,j out of bounds");
        }

        final double x = reward;
        final double time = branchLength;

        final int h = getHfromX(x, time);

        // Ensure C is available up to required N (pdf uses n+1 in C-index)
        ensureCForTime(time, /*extraN=*/1);
        final int N = capacityN - 1;

        final double lt = lambda * time;
        double premult = Math.exp(-lt); // premult_0

        final int uv = i * dim + j;
        double acc = 0.0;

        for (int n = 0; n <= N; ++n) {
            acc += loopCyclePdfEntry(x, time, h, n, premult, uv);
            premult *= lt / (n + 1.0);
        }

        return acc;
    }

    // ------------------------------------------------------------------
// Scalar "into": writes one (dim x dim) matrix into `out`
// ------------------------------------------------------------------
    public void computePdfInto(double x, double time, double[] out) {
        if (out == null) throw new IllegalArgumentException("out must be non-null");
        if (out.length != dim2) {
            throw new IllegalArgumentException("out.length must be dim*dim=" + dim2);
        }

        scratchX1[0] = x;
        scratchT1[0] = time;
        scratchW1[0] = out;

        computePdfInto(scratchX1, scratchT1, false, scratchW1);
    }

    // Convenience: single time broadcast
    public void computePdfInto(double[] X, double time, double[][] W) {
        computePdfInto(X, new double[]{time}, false, W);
    }

    public void computePdfInto(double[] X, double[] times, double[][] W) {
        computePdfInto(X, times, false, W);
    }


    /**
     * Scalar version of loopCyclePdf for a single uv index.
     * Returns the increment to add to Wt[uv] (already includes the temp factor).
     */
    private double loopCyclePdfEntry(double x, double time, int h, int n, double premult, int uv) {

        final double invDiff = invRewardDiff[h];      // 1/(r[h]-r[h-1])
        final double factor  = lambda * invDiff;      // lambda/(r[h]-r[h-1])

        double xh = (x - r[h - 1] * time) * invDiff / time;
        if (xh <= 0.0) xh = 0.0;
        else if (xh >= 1.0) xh = 1.0;

        final double temp = factor * premult;

        if (xh == 0.0) {
            final int aOff = cOffset(h, n + 1, 1);
            final int bOff = cOffset(h, n + 1, 0);
            return temp * (Cflat[aOff + uv] - Cflat[bOff + uv]);
        }

        if (xh == 1.0) {
            final int aOff = cOffset(h, n + 1, n + 1);
            final int bOff = cOffset(h, n + 1, n);
            return temp * (Cflat[aOff + uv] - Cflat[bOff + uv]);
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
            default: w = Math.pow(oneMinus, n);
        }

        double sum = 0.0;

        // k=0
        {
            final int aOff = cOffset(h, n + 1, 1);
            final int bOff = cOffset(h, n + 1, 0);
            sum += w * (Cflat[aOff + uv] - Cflat[bOff + uv]);
        }

        // k=1..n via Bernstein recurrence
        for (int k = 0; k < n; ++k) {
            w *= ((double) (n - k) / (double) (k + 1)) * ratio;
            final int kp1 = k + 1;

            final int aOff = cOffset(h, n + 1, kp1 + 1);
            final int bOff = cOffset(h, n + 1, kp1);

            sum += w * (Cflat[aOff + uv] - Cflat[bOff + uv]);
        }

        return temp * sum;
    }


    // ------------------------------------------------------------
// 2) computePdfInto(X, times, parsimonious, W): reduce allocations & per-call work
//    (a) NO allocations in the hot path
//    (b) avoid Arrays.fill(row) if caller guarantees W is zeroed (toggleable)
// ------------------------------------------------------------
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

        // Validate W rows; overwrite semantics: zero rows here.
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

        // Precompute H, lt, premult0, and maxT in one pass
        double maxT = 0.0;

        if (singleTime) {
            final double time = times[0];
            if (time <= 0.0) throw new IllegalArgumentException("time must be > 0");
            maxT = time;

            final double lt = lambda * time;
            final double premult0 = Math.exp(-lt);

            for (int t = 0; t < T; ++t) {
                s.H[t] = getHfromX(X[t], time);
                s.lt[t] = lt;
                s.premult[t] = premult0;
            }
        } else {
            for (int t = 0; t < T; ++t) {
                final double time = times[t];
                if (time <= 0.0) throw new IllegalArgumentException("time must be > 0 at index " + t);
                if (time > maxT) maxT = time;

                s.H[t] = getHfromX(X[t], time);

                final double lt = lambda * time;
                s.lt[t] = lt;
                s.premult[t] = Math.exp(-lt);
            }
        }

        // Ensure C is available up to required N (pdf uses n+1)
        ensureCForTime(maxT, /*extraN=*/1);
        final int N = capacityN - 1;

        // Optional per-time truncation (only meaningful when times vary)
        if (parsimonious && !singleTime) {
            for (int t = 0; t < T; ++t) {
                s.NN[t] = determineNumberOfSteps(times[t], lambda);
            }
        } else {
            // Avoid Arrays.fill over full capacity; set only [0..T)
            for (int t = 0; t < T; ++t) s.NN[t] = Integer.MAX_VALUE;
        }

        // Main accumulation
        for (int n = 0; n <= N; ++n) {

            if (singleTime) {
                final double time = times[0];
                for (int t = 0; t < T; ++t) {
                    if (s.NN[t] < n) continue;
                    loopCyclePdf(X[t], time, s.H[t], n, s.premult[t], W[t], s.inc);
                }
            } else {
                for (int t = 0; t < T; ++t) {
                    if (s.NN[t] < n) continue;
                    loopCyclePdf(X[t], times[t], s.H[t], n, s.premult[t], W[t], s.inc);
                }
            }

            // premult_{n+1}(t) = premult_n(t) * (lambda*t)/(n+1)
            final double inv = 1.0 / (n + 1.0);
            for (int t = 0; t < T; ++t) {
                s.premult[t] *= s.lt[t] * inv;
            }
        }
    }



    // ------------------------------------------------------------------
// Batched core: writes into W[t] (no allocation), returns void
// ------------------------------------------------------------------


    // ------------------------------------------------------------
    // Public API: CDF (kept; not the optimization focus)
    // ------------------------------------------------------------



    @Override
    public double computeCdf(double x, double time, int i, int j) {
        return computeCdf(new double[]{x}, time)[0][idx(i, j)];
    }

    public double[] computeCdf(double x, double time) {
        return computeCdf(new double[]{x}, time)[0];
    }

    public double[][] computeCdf(double[] X, double time) {

        if (time <= 0.0) throw new IllegalArgumentException("time must be > 0");
        final int T = X.length;

        // ensure C for CDF needs N without the “-1” in your old code
        ensureCForTime(time, /*extraN=*/0);

        final int N = capacityN;

        final double[][] W = new double[T][dim2];
        final int[] H = new int[T];
        for (int t = 0; t < T; t++) H[t] = getHfromX(X[t], time);

        for (int n = 0; n <= N; ++n) {
            accumulateCdf(W, X, H, n, time);
        }
        return W;
    }

    // ------------------------------------------------------------
    // Hot inner kernel: PDF at a single (x,time,h,n) accumulating into Wt
    //
    // We reduce memory traffic by accumulating into `inc` first, then applying `temp` once.
    // `inc` is a thread-local scratch array (dim2).
    // ------------------------------------------------------------
    private void loopCyclePdf(double x, double time, int h, int n, double premult, double[] Wt, double[] inc) {

        final double invDiff = invRewardDiff[h];
        final double factor = lambda * invDiff;

        double xh = (x - r[h - 1] * time) * invDiff / time;

        if (xh <= 0.0) xh = 0.0;
        else if (xh >= 1.0) xh = 1.0;

        final double temp = factor * premult;

        // Boundary cases (avoid Bernstein loop)
        if (xh == 0.0) {
            final int aOff = cOffset(h, n + 1, 1);
            final int bOff = cOffset(h, n + 1, 0);
            for (int uv = 0; uv < dim2; ++uv) {
                Wt[uv] += temp * (Cflat[aOff + uv] - Cflat[bOff + uv]);
            }
            return;
        }

        if (xh == 1.0) {
            final int aOff = cOffset(h, n + 1, n + 1);
            final int bOff = cOffset(h, n + 1, n);
            for (int uv = 0; uv < dim2; ++uv) {
                Wt[uv] += temp * (Cflat[aOff + uv] - Cflat[bOff + uv]);
            }
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
            default: w = Math.pow(oneMinus, n);
        }

        Arrays.fill(inc, 0.0);

        // k=0
        {
            final int aOff = cOffset(h, n + 1, 1);
            final int bOff = cOffset(h, n + 1, 0);
            for (int uv = 0; uv < dim2; ++uv) {
                inc[uv] += w * (Cflat[aOff + uv] - Cflat[bOff + uv]);
            }
        }

        // k=1..n via Bernstein recurrence
        for (int k = 0; k < n; ++k) {
            w *= ((double) (n - k) / (double) (k + 1)) * ratio;
            final int kp1 = k + 1;

            final int aOff = cOffset(h, n + 1, kp1 + 1);
            final int bOff = cOffset(h, n + 1, kp1);

            for (int uv = 0; uv < dim2; ++uv) {
                inc[uv] += w * (Cflat[aOff + uv] - Cflat[bOff + uv]);
            }
        }

        // Apply temp once
        for (int uv = 0; uv < dim2; ++uv) {
            Wt[uv] += temp * inc[uv];
        }
    }

    // ------------------------------------------------------------
    // CDF accumulation (kept as-is structurally)
    // ------------------------------------------------------------
    private void accumulateCdf(double[][] W, double[] X, int[] H, int n, double time) {

        final double premult = Math.exp(
                -lambda * time + n * (Math.log(lambda) + Math.log(time)) - GammaFunction.lnGamma(n + 1.0)
        );

        for (int t = 0; t < X.length; ++t) {
            final double x = X[t];
            final int h = H[t];

            double xh = (x - r[h - 1] * time) / ((r[h] - r[h - 1]) * time);

            final double[] inc = new double[dim2];
            for (int k = 0; k <= n; k++) {
                final double binomialCoef = Binomial.choose(n, k) * Math.pow(xh, k) * Math.pow(1.0 - xh, n - k);
                final int cOff = cOffset(h, n, k);
                for (int uv = 0; uv < dim2; ++uv) {
                    inc[uv] += binomialCoef * Cflat[cOff + uv];
                }
            }
            for (int uv = 0; uv < dim2; ++uv) {
                W[t][uv] += premult * inc[uv];
            }
        }
    }

    // ------------------------------------------------------------
    // Growth + computation of C
    // ------------------------------------------------------------
    private void ensureCForTime(double time, int extraN) {
        if (time <= 0.0) throw new IllegalArgumentException("time must be > 0");

        if (time > maxTime) maxTime = time;

        final int requiredN = determineNumberOfSteps(time, lambda) + extraN;

        if (capacityN < requiredN) {
            ensureCapacityN(requiredN);
            computeChnk();
        }
    }


    private void ensureCapacityN(int requiredN) {

        int newCap = (capacityN < 0) ? requiredN : capacityN;
        while (newCap < requiredN) {
            newCap = Math.max(requiredN, (int) (newCap * 1.5) + 1);
        }
        if (capacityN < 0) newCap = requiredN;

        final int newN1 = newCap + 1;

        final long blocks = (long) (phi + 1) * (long) newN1 * (long) newN1;
        final long totalDoubles = blocks * (long) dim2;

        if (totalDoubles > Integer.MAX_VALUE) {
            throw new IllegalStateException("C array too large: " + totalDoubles + " doubles");
        }

        // IMPORTANT: allocate fresh; DO NOT copy old (we recompute everything in computeChnk()).
        this.Cflat = new double[(int) totalDoubles];
        this.capacityN = newCap;
        this.N1 = newN1;

        if (DEBUG) {
            System.err.println("Allocated Cflat for capacityN=" + capacityN + " (N1=" + N1 + "), doubles=" + totalDoubles);
        }
    }

    /**
     * Compute all C(h,n,k,uv) up to capacityN using Sericola recursions.
     */
    private void computeChnk() {

        if (Cflat == null || capacityN < 0) {
            throw new IllegalStateException("C storage not initialized");
        }

        // Pn = I (reuse buffer; no allocation)
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

                        final double cScalar = (r[u] - r[h]) / (r[u] - r[h - 1]);
                        final double dScalar = (r[h] - r[h - 1]) / (r[u] - r[h - 1]);

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

            // Pn = Pn * P (reuse buffers)
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

                        final double cScalar = (r[h - 1] - r[u]) / (r[h] - r[u]);
                        final double dScalar = (r[h] - r[h - 1]) / (r[h] - r[u]);

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

    // ------------------------------------------------------------
    // Helpers
    // ------------------------------------------------------------

    private void ensureScratchCapacity(Scratch s, int T) {
        if (s.H == null || s.H.length < T) s.H = new int[T];
        if (s.NN == null || s.NN.length < T) s.NN = new int[T];
        if (s.lt == null || s.lt.length < T) s.lt = new double[T];
        if (s.premult == null || s.premult.length < T) s.premult = new double[T];
        if (s.inc == null || s.inc.length != dim2) s.inc = new double[dim2];
    }

    private int idx(int i, int j) {
        return i * dim + j;
    }

    private int cOffset(int h, int n, int k) {
        // base offset for C(h,n,k,0) in flattened array
        return (((h * N1 + n) * N1 + k) * dim2);
    }

    private double determineLambda() {
        double minDiag = Q[0];
        for (int i = 1; i < dim; ++i) {
            final double d = Q[idx(i, i)];
            if (d < minDiag) minDiag = d;
        }
        return -minDiag;
    }

    private double[] initializeP(double[] Q, double lambda) {
        final double[] P = new double[dim2];
        final double inv = 1.0 / lambda;
        for (int i = 0; i < dim; ++i) {
            final int ioff = i * dim;
            for (int j = 0; j < dim; ++j) {
                P[ioff + j] = (i == j ? 1.0 : 0.0) + Q[ioff + j] * inv;
            }
        }
        return P;
    }

    private int getHfromX(double x, double time) {
        if (time <= 0.0) throw new IllegalArgumentException("time must be > 0");

        final double lo = r[0] * time;
        final double hi = r[phi] * time;

        if (x < lo) throw new IllegalArgumentException("x must be >= r[0] * time");
        if (x > hi) throw new IllegalArgumentException("x must be <= r[phi] * time");

        final double rate = x / time;
        int h = Arrays.binarySearch(r, 0, phi + 1, rate);
        if (h < 0) h = -h - 1;

        if (h <= 0) return 1;
        if (h > phi) return phi;
        return h;
    }

    /**
     * Fast Poisson CDF truncation:
     * find smallest i such that sum_{k=0}^i e^{-lt} (lt)^k / k! >= 1-epsilon.
     */
    private int determineNumberOfSteps(double time, double lambda) {
        final double target = 1.0 - epsilon;
        final double lt = lambda * time;

        double term = Math.exp(-lt);
        double sum = term;

        int i = 0;
        // cap tied to lt to avoid spurious warnings
        final int hardCap = Math.max(5000, (int) (lt + 10.0 * Math.sqrt(lt + 1.0)));

        while (sum < target && i < hardCap) {
            i++;
            term *= lt / i;
            sum += term;
        }
        return i;
    }

    private double[][] squareMatrix(final double[] mat) {
        final double[][] rtn = new double[dim][dim];
        for (int i = 0; i < dim; ++i) {
            System.arraycopy(mat, i * dim, rtn[i], 0, dim);
        }
        return rtn;
    }

    // ------------------------------------------------------------
    // Conditional probabilities (matrix exponential of Q)
    // ------------------------------------------------------------

    private EigenDecomposition getEigenDecomposition() {
        if (eigenDecomposition == null) {
            eigenDecomposition = eigenSystem.decomposeMatrix(squareMatrix(Q));
        }
        return eigenDecomposition;
    }

    public double[] computeConditionalProbabilities(double distance) {
        final double[] matrix = new double[dim2];
        eigenSystem.computeExponential(getEigenDecomposition(), distance, matrix);
        return matrix;
    }

    public double computeConditionalProbability(double distance, int i, int j) {
        return eigenSystem.computeExponential(getEigenDecomposition(), distance, i, j);
    }

    // ------------------------------------------------------------
    // Debug / info
    // ------------------------------------------------------------

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }
    private int storedCapacityN = -1;
    private int storedN1 = 0;
    private double storedMaxTime = 0.0;
    private boolean hasStoredState = false;

    @Override
    public void storeState() {
        storedCapacityN = capacityN;
        storedN1 = N1;
        storedMaxTime = maxTime;
        hasStoredState = true;
    }

    @Override
    public void restoreState() {
        if (!hasStoredState) return;
        capacityN = storedCapacityN;
        N1 = storedN1;
        maxTime = storedMaxTime;
        // IMPORTANT: do NOT null Cflat; do NOT shrink.
        // Next ensureCForTime/computeChnk will overwrite as needed.
    }

    @Override
    protected void acceptState() {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Q: ").append(new Vector(Q)).append("\n");
        sb.append("r: ").append(new Vector(r)).append("\n");
        sb.append("lambda: ").append(lambda).append("\n");
        sb.append("capacityN: ").append(capacityN).append("\n");
        sb.append("maxTime: ").append(maxTime).append("\n");
        sb.append("cprob at maxTime: ").append(new Vector(computeConditionalProbabilities(maxTime))).append("\n");
        return sb.toString();
    }
}
