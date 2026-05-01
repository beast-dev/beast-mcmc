package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalNumerics;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import org.ejml.data.DenseMatrix64F;

/**
 * Branch-local canonical transition construction for OU tree gradients.
 */
final class OUCanonicalTransitionBuilder {

    private OUCanonicalTransitionBuilder() { }

    static void fillFromPrecisionMoments(final DenseMatrix64F actualization,
                                         final DenseMatrix64F displacement,
                                         final DenseMatrix64F precision,
                                         final CanonicalGaussianTransition out) {
        final int dimension = out.getDimension();
        for (int i = 0; i < dimension; ++i) {
            double infoY = 0.0;
            final int iOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                final double p = precision.unsafe_get(i, j);
                out.precisionYY[iOffset + j] = p;
                out.precisionYX[iOffset + j] = -multiplyEntry(precision, actualization, i, j);
                out.precisionXY[j * dimension + i] = out.precisionYX[iOffset + j];
                out.precisionXX[iOffset + j] =
                        multiplyEntryTranspose(actualization, precision, actualization, i, j);
                infoY += p * displacement.unsafe_get(j, 0);
            }
            out.informationY[i] = infoY;
        }
        for (int i = 0; i < dimension; ++i) {
            double infoX = 0.0;
            for (int j = 0; j < dimension; ++j) {
                infoX -= actualization.unsafe_get(j, i) * out.informationY[j];
            }
            out.informationX[i] = infoX;
        }
        out.logNormalizer = 0.0;
    }

    static void fillFromMoments(final double[] transitionMatrix,
                                final double[] transitionOffset,
                                final double[] transitionCovariance,
                                final CanonicalGaussianTransition out,
                                final double[] precisionScratch,
                                final double[] symmetricScratch,
                                final double[] choleskyScratch,
                                final double[] lowerInverseScratch,
                                final double[] transitionTransposeScratch,
                                final double[] matrixScratch) {
        final int dimension = out.getDimension();
        final double logDet = invertSymmetricPositiveDefinite(
                transitionCovariance,
                dimension,
                precisionScratch,
                symmetricScratch,
                choleskyScratch,
                lowerInverseScratch);
        MatrixOps.transpose(transitionMatrix, transitionTransposeScratch, dimension);

        MatrixOps.matMul(transitionTransposeScratch, precisionScratch, matrixScratch, dimension);
        MatrixOps.matMul(matrixScratch, transitionMatrix, out.precisionXX, dimension);

        MatrixOps.matMul(transitionTransposeScratch, precisionScratch, out.precisionXY, dimension);
        MatrixOps.scaleInPlace(out.precisionXY, -1.0, dimension * dimension);

        MatrixOps.matMul(precisionScratch, transitionMatrix, out.precisionYX, dimension);
        MatrixOps.scaleInPlace(out.precisionYX, -1.0, dimension * dimension);

        System.arraycopy(precisionScratch, 0, out.precisionYY, 0, dimension * dimension);

        MatrixOps.matVec(precisionScratch, transitionOffset, out.informationY, dimension);
        MatrixOps.matVec(transitionTransposeScratch, out.informationY, out.informationX, dimension);
        MatrixOps.scaleInPlace(out.informationX, -1.0, dimension);

        out.logNormalizer =
                0.5 * (dimension * Math.log(2.0 * Math.PI)
                        + logDet
                        + dot(transitionOffset, out.informationY, dimension));
    }

    private static double multiplyEntry(final DenseMatrix64F left,
                                        final DenseMatrix64F right,
                                        final int row,
                                        final int column) {
        double sum = 0.0;
        final int inner = left.numCols;
        for (int k = 0; k < inner; ++k) {
            sum += left.unsafe_get(row, k) * right.unsafe_get(k, column);
        }
        return sum;
    }

    private static double multiplyEntryTranspose(final DenseMatrix64F left,
                                                 final DenseMatrix64F middle,
                                                 final DenseMatrix64F right,
                                                 final int row,
                                                 final int column) {
        double sum = 0.0;
        final int inner = left.numRows;
        for (int i = 0; i < inner; ++i) {
            double leftMiddle = 0.0;
            for (int k = 0; k < inner; ++k) {
                leftMiddle += left.unsafe_get(k, row) * middle.unsafe_get(k, i);
            }
            sum += leftMiddle * right.unsafe_get(i, column);
        }
        return sum;
    }

    private static double dot(final double[] left, final double[] right, final int dimension) {
        double sum = 0.0;
        for (int i = 0; i < dimension; ++i) {
            sum += left[i] * right[i];
        }
        return sum;
    }

    private static double invertSymmetricPositiveDefinite(final double[] matrix,
                                                          final int dimension,
                                                          final double[] inverseOut,
                                                          final double[] symmetricScratch,
                                                          final double[] choleskyScratch,
                                                          final double[] lowerInverseScratch) {
        final int d2 = dimension * dimension;
        for (int i = 0; i < dimension; ++i) {
            final int iOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                symmetricScratch[iOffset + j] = 0.5 * (matrix[iOffset + j] + matrix[j * dimension + i]);
            }
        }

        final double jitterBase = Math.max(
                1.0e-12,
                1.0e-12 * Math.max(1.0, maxAbsDiagonal(symmetricScratch, dimension)));
        double jitter = 0.0;
        for (int attempt = 0; attempt < 12; ++attempt) {
            System.arraycopy(symmetricScratch, 0, choleskyScratch, 0, d2);
            if (jitter > 0.0) {
                for (int i = 0; i < dimension; ++i) {
                    choleskyScratch[i * dimension + i] += jitter;
                }
            }
            if (MatrixOps.tryCholesky(choleskyScratch, lowerInverseScratch, dimension)) {
                return MatrixOps.invertFromCholesky(lowerInverseScratch, choleskyScratch, inverseOut, dimension);
            }
            jitter = jitter == 0.0 ? jitterBase : jitter * 10.0;
        }
        throw new IllegalArgumentException("Matrix is not symmetric positive definite.");
    }

    private static double maxAbsDiagonal(final double[] matrix, final int dimension) {
        double max = 0.0;
        for (int i = 0; i < dimension; ++i) {
            max = Math.max(max, Math.abs(matrix[i * dimension + i]));
        }
        return max;
    }
}
