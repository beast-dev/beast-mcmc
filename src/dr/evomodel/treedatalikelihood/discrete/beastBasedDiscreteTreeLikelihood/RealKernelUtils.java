package dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood;

import dr.evomodel.substmodel.EigenDecomposition;
import org.apache.commons.math.util.FastMath;

/**
 * Fréchet-integral kernel for substitution models with all-real eigenvalues
 * (e.g. reversible CTMCs satisfying detailed balance).
 *
 * For an all-real eigensystem every block is 1×1 (scalar), so the Fréchet
 * integral coefficient for entry (i, j) is:
 *
 *   c_{ij}(t) = (exp(t λᵢ) − exp(t λⱼ)) / (λᵢ − λⱼ)   λᵢ ≠ λⱼ
 *             = t · exp(t λᵢ)                             λᵢ = λⱼ
 *
 * Compared to {@link ComplexBlockKernelUtils} this class:
 *   - Eliminates all cos / sin calls.
 *   - Uses a packed K² coefficient array (stride 1) instead of K²×16
 *     (stride 16).  For K = 17 that shrinks the hot working set from ~37 KB
 *     (L2/L3) to ~2.3 KB (L1-resident).
 *   - Removes all block-type dispatch overhead.
 */
public final class RealKernelUtils {

    private static final double REAL_TOLERANCE = 1.0e-12;

    private RealKernelUtils() {}

    // -------------------------------------------------------------------------
    // Plan — all storage pre-allocated at construction time, zero heap allocation
    // during gradient evaluation.
    // -------------------------------------------------------------------------

    public static final class RealKernelPlan {

        /** Time-independent: 1/(λᵢ − λⱼ), or 0 when |λᵢ − λⱼ| < tol. Layout: [i*K+j]. */
        final double[] invDenom;

        /** Time-dependent: exp(t·λᵢ). Layout: [i]. */
        final double[] expR;

        /**
         * Time-dependent packed coefficients c_{ij}(t). Layout: [i*K+j].
         * Only populated by {@link #fillCoefficients} (the multi-pattern path).
         */
        final double[] coefficients;

        public RealKernelPlan(int stateCount) {
            this.invDenom     = new double[stateCount * stateCount];
            this.expR         = new double[stateCount];
            this.coefficients = new double[stateCount * stateCount];
        }
    }

    // -------------------------------------------------------------------------
    // Detection
    // -------------------------------------------------------------------------

    /**
     * Returns true when every imaginary part of the eigenvalues is below
     * {@code REAL_TOLERANCE}, i.e. the model has an all-real eigensystem.
     */
    public static boolean isAllReal(EigenDecomposition eigenDecomp, int stateCount) {
        final double[] ev = eigenDecomp.getEigenValues();
        if (ev.length == stateCount) return true;
        for (int i = 0; i < stateCount; i++) {
            if (Math.abs(ev[stateCount + i]) > REAL_TOLERANCE) return false;
        }
        return true;
    }

    // -------------------------------------------------------------------------
    // Time-independent structure — called once per eigensystem update.
    // -------------------------------------------------------------------------

    /**
     * Fills the time-independent inverse-denominator table.
     * Must be called whenever the substitution model parameters change.
     */
    public static void fillStructure(RealKernelPlan plan,
                                      EigenDecomposition eigenDecomp,
                                      int stateCount) {
        final double[] ev = eigenDecomp.getEigenValues();
        for (int i = 0; i < stateCount; i++) {
            final double li       = ev[i];
            final int    rowBase  = i * stateCount;
            for (int j = 0; j < stateCount; j++) {
                final double delta = li - ev[j];
                plan.invDenom[rowBase + j] =
                        Math.abs(delta) < REAL_TOLERANCE ? 0.0 : 1.0 / delta;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Single-pattern fused fill+apply (patternCount == 1, e.g. phylogeography).
    // Computes coefficients on the fly — no K²-element intermediate store needed.
    // -------------------------------------------------------------------------

    public static void fillAndApplyToOuterProduct(RealKernelPlan plan,
                                                   EigenDecomposition eigenDecomp,
                                                   double time,
                                                   double[] leftVector,
                                                   double[] rightVector,
                                                   double scale,
                                                   double[] eigenBasisGradient,
                                                   int stateCount) {
        final double[] ev      = eigenDecomp.getEigenValues();
        final double[] expR    = plan.expR;
        final double[] invD    = plan.invDenom;

        for (int i = 0; i < stateCount; i++) {
            expR[i] = FastMath.exp(time * ev[i]);
        }

        for (int i = 0; i < stateCount; i++) {
            final double li = scale * leftVector[i];
            if (li == 0.0) continue;
            final double expI   = expR[i];
            final int    rowBase = i * stateCount;
            for (int j = 0; j < stateCount; j++) {
                final double d     = invD[rowBase + j];
                final double coeff = d == 0.0 ? time * expI : (expI - expR[j]) * d;
                eigenBasisGradient[rowBase + j] += coeff * li * rightVector[j];
            }
        }
    }

    // -------------------------------------------------------------------------
    // Multi-pattern path: fill coefficients once, apply across all patterns.
    // -------------------------------------------------------------------------

    public static void fillCoefficients(RealKernelPlan plan,
                                         EigenDecomposition eigenDecomp,
                                         double time,
                                         int stateCount) {
        final double[] ev   = eigenDecomp.getEigenValues();
        final double[] expR = plan.expR;
        final double[] invD = plan.invDenom;
        final double[] c    = plan.coefficients;

        for (int i = 0; i < stateCount; i++) {
            expR[i] = FastMath.exp(time * ev[i]);
        }
        for (int i = 0; i < stateCount; i++) {
            final double expI   = expR[i];
            final int    rowBase = i * stateCount;
            for (int j = 0; j < stateCount; j++) {
                final double d = invD[rowBase + j];
                c[rowBase + j] = d == 0.0 ? time * expI : (expI - expR[j]) * d;
            }
        }
    }

    public static void applyToOuterProduct(RealKernelPlan plan,
                                            double[] leftVector,
                                            double[] rightVector,
                                            double scale,
                                            double[] eigenBasisGradient,
                                            int stateCount) {
        final double[] c = plan.coefficients;
        for (int i = 0; i < stateCount; i++) {
            final double li = scale * leftVector[i];
            if (li == 0.0) continue;
            final int rowBase = i * stateCount;
            for (int j = 0; j < stateCount; j++) {
                eigenBasisGradient[rowBase + j] += c[rowBase + j] * li * rightVector[j];
            }
        }
    }
}
