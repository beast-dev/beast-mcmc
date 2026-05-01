package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalNumerics;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalNumericsOptions;
import org.ejml.data.DenseMatrix64F;

/**
 * Positive-definite inversion helper used by orthogonal block canonical OU code.
 */
final class OrthogonalBlockPositiveDefiniteInverter {

    private static final CanonicalNumericsOptions OPTIONS = CanonicalNumericsOptions.ORTHOGONAL_BLOCK;

    private OrthogonalBlockPositiveDefiniteInverter() { }

    static double copyAndInvert(final DenseMatrix64F source,
                                final double[][] matrixOut,
                                final double[][] inverseOut,
                                final double[] choleskyScratch,
                                final double[] lowerInverseScratch) {
        return CanonicalNumerics.copyAndInvertPositiveDefinite(
                source, matrixOut, inverseOut, choleskyScratch, lowerInverseScratch, OPTIONS);
    }

    static double copyAndInvertFlat(final DenseMatrix64F source,
                                    final double[] matrixOut,
                                    final double[] inverseOut,
                                    final double[] choleskyScratch,
                                    final double[] lowerInverseScratch) {
        return CanonicalNumerics.copyAndInvertPositiveDefiniteFlat(
                source, matrixOut, inverseOut, choleskyScratch, lowerInverseScratch, OPTIONS);
    }
}
