package dr.evomodel.continuous.ou.blockdiagonal;

import dr.evomodel.continuous.ou.math.OUPositiveDefiniteMatrixInversionOptions;
import dr.evomodel.continuous.ou.math.PositiveDefiniteMatrixOps;
import org.ejml.data.DenseMatrix64F;

/**
 * Positive-definite inversion helper used by orthogonal block canonical OU code.
 */
final class BlockDiagonalPositiveDefiniteInverter {

    private static final OUPositiveDefiniteMatrixInversionOptions OPTIONS =
            OUPositiveDefiniteMatrixInversionOptions.ORTHOGONAL_BLOCK;

    private BlockDiagonalPositiveDefiniteInverter() { }

    static double copyAndInvertFlat(final DenseMatrix64F source,
                                    final double[] matrixOut,
                                    final double[] inverseOut,
                                    final double[] choleskyScratch,
                                    final double[] lowerInverseScratch) {
        return PositiveDefiniteMatrixOps.copyAndInvertPositiveDefiniteFlat(
                source, matrixOut, inverseOut, choleskyScratch, lowerInverseScratch, OPTIONS);
    }
}
