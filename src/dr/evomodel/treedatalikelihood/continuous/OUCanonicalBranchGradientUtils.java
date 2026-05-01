package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
final class OUCanonicalBranchGradientUtils {

    private OUCanonicalBranchGradientUtils() {
        // no instances
    }

    static boolean isFinite(final double[] values) {
        for (double value : values) {
            if (!Double.isFinite(value)) {
                return false;
            }
        }
        return true;
    }

    static void zero(final double[] vector) {
        for (int i = 0; i < vector.length; ++i) {
            vector[i] = 0.0;
        }
    }

    static void transposeFlatSquare(final double[] in, final double[] out, final int dim) {
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                out[j * dim + i] = in[i * dim + j];
            }
        }
    }

    static void copyAdjoints(final CanonicalLocalTransitionAdjoints source,
                             final CanonicalLocalTransitionAdjoints target,
                             final boolean includeTransition,
                             final boolean includeCovariance) {
        final int matLen = source.dLogL_dF.length;
        final int vecLen = source.dLogL_df.length;
        for (int k = 0; k < matLen; ++k) {
            target.dLogL_dF[k] = includeTransition ? source.dLogL_dF[k] : 0.0;
            target.dLogL_dOmega[k] = includeCovariance ? source.dLogL_dOmega[k] : 0.0;
        }
        for (int i = 0; i < vecLen; ++i) {
            target.dLogL_df[i] = includeTransition ? source.dLogL_df[i] : 0.0;
        }
    }
}
