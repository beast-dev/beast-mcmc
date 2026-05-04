package dr.evomodel.treedatalikelihood.continuous.canonical.math;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;

/**
 * Converts between moment-form and canonical-form Gaussian representations.
 *
 * <p>Primary methods use flat row-major arrays and caller-owned {@link Workspace}
 * instances. The {@code double[][]} overloads are boundary helpers for BEAST model
 * objects and allocate scratch.
 *
 * <h3>Numerical robustness</h3>
 * <p>Inversion uses Cholesky factorization. If the matrix is not strictly positive
 * definite, up to {@value #MAX_JITTER_ATTEMPTS} retries are made with increasing
 * diagonal jitter. A final failure throws {@link IllegalArgumentException}.
 */
public final class GaussianFormConverter {

    private static final double JITTER_RELATIVE = 1.0e-12;
    private static final double JITTER_ABSOLUTE  = 1.0e-12;
    private static final int    MAX_JITTER_ATTEMPTS = 12;

    private GaussianFormConverter() { }

    // -----------------------------------------------------------------------
    // Workspace
    // -----------------------------------------------------------------------

    /**
     * Reusable scratch buffers for conversion methods.
     *
     * <p>Call {@link #ensureDim(int)} before passing to a conversion method.
     * A single {@code Workspace} may be shared across repeated calls to any method
     * in this class, provided the calls are not concurrent.
     */
    public static final class Workspace {

        public double[] scratch    = new double[0]; // working copy of input, or general temp matrix
        public double[] lower      = new double[0]; // Cholesky lower factor L, or P = Ω⁻¹
        public double[] lowerInv   = new double[0]; // L⁻¹, or F^T (reused after inversion)
        public double[] temp       = new double[0]; // general temp matrix (e.g. F^T P)
        public double[] symmetric  = new double[0]; // symmetrized copy for jitter retry
        public double[] adjusted   = new double[0]; // symmetrized + diagonal jitter
        public double[] tempVec    = new double[0]; // dim-length scratch vector
        private int dim            = 0;

        public Workspace() { }

        /**
         * Ensures all scratch arrays are sized for matrices of dimension {@code newDim}.
         * No-op if already large enough.
         */
        public void ensureDim(int newDim) {
            if (newDim <= dim) return;
            final int d2 = newDim * newDim;
            scratch    = new double[d2];
            lower      = new double[d2];
            lowerInv   = new double[d2];
            temp       = new double[d2];
            symmetric  = new double[d2];
            adjusted   = new double[d2];
            tempVec    = new double[newDim];
            dim        = newDim;
        }
    }

    // -----------------------------------------------------------------------
    // State: moments → canonical
    // -----------------------------------------------------------------------

    /**
     * Fills {@code out} with the canonical-form Gaussian state for N(μ, Σ).
     *
     * <p>Sets:
     * <pre>
     *   J = Σ⁻¹,   h = J μ,   g = ½(d log 2π + log|Σ| + μ^T h)
     * </pre>
     *
     * @param mean    mean vector μ, length {@code dim}
     * @param flatCov row-major covariance Σ, length {@code dim*dim}
     * @param dim     state dimension
     * @param ws      caller-owned workspace; call {@code ws.ensureDim(dim)} first
     * @param out     canonical state to fill; must have dimension {@code dim}
     */
    public static void fillStateFromMoments(
            double[] mean, double[] flatCov, int dim, Workspace ws, CanonicalGaussianState out) {
        System.arraycopy(flatCov, 0, ws.scratch, 0, dim * dim);
        final double logDet = invertSPDRobust(ws.scratch, out.precision, dim, ws);
        MatrixOps.matVec(out.precision, mean, out.information, dim);
        out.logNormalizer = 0.5 * (dim * MatrixOps.LOG_TWO_PI
                + logDet + dot(mean, out.information, dim));
    }

    /**
     * Boundary overload accepting a {@code double[][]} covariance.
     * Allocates a {@code Workspace}; prefer {@link #fillStateFromMoments(double[], double[], int, Workspace, CanonicalGaussianState)}
     * on hot paths.
     */
    public static void fillStateFromMoments(
            double[] mean, double[][] covariance, CanonicalGaussianState out) {
        final int dim = mean.length;
        final Workspace ws = new Workspace();
        ws.ensureDim(dim);
        MatrixOps.toFlat(covariance, ws.scratch, dim);
        fillStateFromMoments(mean, ws.scratch, dim, ws, out);
    }

    // -----------------------------------------------------------------------
    // Transition: moments → canonical
    // -----------------------------------------------------------------------

    /**
     * Fills {@code out} with the canonical-form Gaussian transition for
     * p(y | x) = N(F x + f, Ω).
     *
     * <p>The joint block-precision over [x, y] is:
     * <pre>
     *   J_xx =  F^T Ω⁻¹ F,   J_xy = J_yx^T = -F^T Ω⁻¹
     *   J_yy =  Ω⁻¹
     *   h_x  = -F^T Ω⁻¹ f,   h_y  = Ω⁻¹ f
     *   g    =  ½(d log 2π + log|Ω| + f^T Ω⁻¹ f)
     * </pre>
     *
     * @param flatF    row-major transition matrix F, length {@code dim*dim};
     *                 must not alias {@code ws.lowerInv} or {@code ws.temp}
     * @param f        transition offset vector, length {@code dim}
     * @param flatOmega row-major transition covariance Ω, length {@code dim*dim}
     * @param dim      state dimension
     * @param ws       caller-owned workspace; call {@code ws.ensureDim(dim)} first
     * @param out      canonical transition to fill; must have dimension {@code dim}
     */
    public static void fillTransitionFromMoments(
            double[] flatF, double[] f, double[] flatOmega,
            int dim, Workspace ws, CanonicalGaussianTransition out) {

        // Invert Ω → P = Ω⁻¹, stored in ws.lower (scratch holds Ω, untouched by Cholesky)
        System.arraycopy(flatOmega, 0, ws.scratch, 0, dim * dim);
        final double logDet = invertSPDRobust(ws.scratch, ws.lower, dim, ws);
        // ws.lower = P = Ω⁻¹

        // F^T → ws.lowerInv
        MatrixOps.transpose(flatF, ws.lowerInv, dim);

        // F^T P → ws.temp
        MatrixOps.matMul(ws.lowerInv, ws.lower, ws.temp, dim);

        // precisionXX = F^T P F
        MatrixOps.matMul(ws.temp, flatF, out.precisionXX, dim);

        // precisionXY = -F^T P  (copy ws.temp, negate)
        System.arraycopy(ws.temp, 0, out.precisionXY, 0, dim * dim);
        MatrixOps.scaleInPlace(out.precisionXY, -1.0, dim * dim);

        // precisionYX = -P F = -(F^T P)^T  (transpose ws.temp, negate)
        MatrixOps.transpose(ws.temp, out.precisionYX, dim);
        MatrixOps.scaleInPlace(out.precisionYX, -1.0, dim * dim);

        // precisionYY = P
        System.arraycopy(ws.lower, 0, out.precisionYY, 0, dim * dim);

        // informationY = P f
        MatrixOps.matVec(ws.lower, f, out.informationY, dim);

        // informationX = -F^T (P f) = -F^T * informationY
        MatrixOps.matVec(ws.lowerInv, out.informationY, out.informationX, dim);
        MatrixOps.scaleInPlace(out.informationX, -1.0, dim);

        out.logNormalizer = 0.5 * (dim * MatrixOps.LOG_TWO_PI
                + logDet + dot(f, out.informationY, dim));
    }

    /**
     * Boundary overload accepting {@code double[][]} F and Omega.
     * Allocates scratch arrays; prefer the flat variant on hot paths.
     */
    public static void fillTransitionFromMoments(
            double[][] F, double[] f, double[][] Omega, CanonicalGaussianTransition out) {
        final int dim = f.length;
        final Workspace ws = new Workspace();
        ws.ensureDim(dim);
        final double[] flatF = new double[dim * dim];
        MatrixOps.toFlat(F, flatF, dim);
        MatrixOps.toFlat(Omega, ws.scratch, dim);
        fillTransitionFromMoments(flatF, f, ws.scratch, dim, ws, out);
    }

    // -----------------------------------------------------------------------
    // State: canonical → moments
    // -----------------------------------------------------------------------

    /**
     * Fills {@code meanOut} and {@code flatCovOut} with the moment-form (μ, Σ)
     * of the canonical-form state {@code state}.
     *
     * <p>Sets:
     * <pre>
     *   Σ = J⁻¹,   μ = Σ h
     * </pre>
     *
     * @param state      canonical Gaussian state
     * @param meanOut    output mean μ, length {@code dim}
     * @param flatCovOut output covariance Σ row-major, length {@code dim*dim}
     * @param dim        state dimension
     * @param ws         caller-owned workspace; call {@code ws.ensureDim(dim)} first
     */
    public static void fillMomentsFromState(
            CanonicalGaussianState state,
            double[] meanOut, double[] flatCovOut, int dim, Workspace ws) {
        System.arraycopy(state.precision, 0, ws.scratch, 0, dim * dim);
        invertSPDRobustNoLogDet(ws.scratch, flatCovOut, dim, ws);
        MatrixOps.matVec(flatCovOut, state.information, meanOut, dim);
    }

    /**
     * Bridge overload writing covariance into a {@code double[][]}.
     * Allocates scratch arrays; prefer the flat variant on hot paths.
     */
    public static void fillMomentsFromState(
            CanonicalGaussianState state, double[] meanOut, double[][] covOut) {
        final int dim = meanOut.length;
        final Workspace ws = new Workspace();
        ws.ensureDim(dim);
        final double[] flatCov = new double[dim * dim];
        fillMomentsFromState(state, meanOut, flatCov, dim, ws);
        MatrixOps.fromFlat(flatCov, covOut, dim);
    }

    // -----------------------------------------------------------------------
    // Private: robust SPD inversion
    // -----------------------------------------------------------------------

    /**
     * Invert a symmetric positive-definite matrix {@code src} into {@code dst},
     * retrying with increasing diagonal jitter if Cholesky factorization fails.
     *
     * <p>{@code src} is not modified. Uses {@code ws.lower}, {@code ws.lowerInv},
     * {@code ws.symmetric}, and {@code ws.adjusted} as scratch.
     *
     * @return log-determinant of {@code src} (with any applied jitter)
     * @throws IllegalArgumentException after all jitter attempts are exhausted
     */
    static double invertSPDRobust(double[] src, double[] dst, int dim, Workspace ws) {
        // Fast path: try strict Cholesky directly on src (src not modified by tryCholesky)
        if (MatrixOps.tryCholesky(src, ws.lower, dim)) {
            return MatrixOps.invertFromCholesky(ws.lower, ws.lowerInv, dst, dim);
        }

        // Symmetrize src → ws.symmetric
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                ws.symmetric[i * dim + j] = 0.5 * (src[i * dim + j] + src[j * dim + i]);

        // Compute Gershgorin lower bound to choose initial jitter
        double maxDiag = 0.0;
        for (int i = 0; i < dim; i++)
            maxDiag = Math.max(maxDiag, Math.abs(ws.symmetric[i * dim + i]));
        final double jitterBase = Math.max(JITTER_ABSOLUTE, JITTER_RELATIVE * Math.max(1.0, maxDiag));

        double jitter = 0.0;
        double lowerBound = Double.NaN;

        for (int attempt = 0; attempt < MAX_JITTER_ATTEMPTS; attempt++) {
            System.arraycopy(ws.symmetric, 0, ws.adjusted, 0, dim * dim);
            if (jitter > 0.0) {
                for (int i = 0; i < dim; i++) ws.adjusted[i * dim + i] += jitter;
            }
            if (MatrixOps.tryCholesky(ws.adjusted, ws.lower, dim)) {
                return MatrixOps.invertFromCholesky(ws.lower, ws.lowerInv, dst, dim);
            }
            if (jitter == 0.0) {
                if (Double.isNaN(lowerBound)) lowerBound = gershgorinLowerBound(ws.symmetric, dim);
                jitter = lowerBound > 0.0 ? jitterBase : (-lowerBound + jitterBase);
            } else {
                jitter *= 10.0;
            }
        }
        throw new IllegalArgumentException(
                "Matrix is not positive definite after " + MAX_JITTER_ATTEMPTS + " jitter retries"
                + "; dim=" + dim + "; maxDiag=" + maxDiag);
    }

    /** Same as {@link #invertSPDRobust} but discards the log-determinant. */
    static void invertSPDRobustNoLogDet(double[] src, double[] dst, int dim, Workspace ws) {
        invertSPDRobust(src, dst, dim, ws);
    }

    // -----------------------------------------------------------------------
    // Private utilities
    // -----------------------------------------------------------------------

    private static double dot(double[] a, double[] b, int dim) {
        double s = 0.0;
        for (int i = 0; i < dim; i++) s += a[i] * b[i];
        return s;
    }

    private static double gershgorinLowerBound(double[] m, int dim) {
        double lower = Double.POSITIVE_INFINITY;
        for (int i = 0; i < dim; i++) {
            double radius = 0.0;
            for (int j = 0; j < dim; j++) if (j != i) radius += Math.abs(m[i * dim + j]);
            lower = Math.min(lower, m[i * dim + i] - radius);
        }
        return lower;
    }
}
