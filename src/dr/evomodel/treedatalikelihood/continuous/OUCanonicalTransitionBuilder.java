package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalNumerics;
import org.ejml.data.DenseMatrix64F;

/**
 * Branch-local transition construction for the time-series/tree canonical bridge.
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

    static void fillFromMoments(final double[][] transitionMatrix,
                                final double[] transitionOffset,
                                final double[][] transitionCovariance,
                                final CanonicalGaussianTransition out,
                                final double[][] precisionScratch,
                                final double[][] symmetricScratch,
                                final double[][] choleskyScratch,
                                final double[][] lowerInverseScratch,
                                final double[][] transitionTransposeScratch,
                                final double[][] matrixScratch) {
        final int dimension = out.getDimension();
        final double logDet = CanonicalNumerics.invertSymmetricPositiveDefinite(
                transitionCovariance,
                dimension,
                precisionScratch,
                symmetricScratch,
                choleskyScratch,
                lowerInverseScratch);
        transpose(transitionMatrix, transitionTransposeScratch, dimension);

        multiply(transitionTransposeScratch, precisionScratch, matrixScratch, dimension);
        multiply(matrixScratch, transitionMatrix, out.precisionXX, dimension);

        multiply(transitionTransposeScratch, precisionScratch, out.precisionXY, dimension);
        scaleInPlace(out.precisionXY, -1.0, dimension * dimension);

        multiply(precisionScratch, transitionMatrix, out.precisionYX, dimension);
        scaleInPlace(out.precisionYX, -1.0, dimension * dimension);

        copyMatrix(precisionScratch, out.precisionYY, dimension);

        multiply(precisionScratch, transitionOffset, out.informationY, dimension);
        multiply(transitionTransposeScratch, out.informationY, out.informationX, dimension);
        scaleInPlace(out.informationX, -1.0, dimension);

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

    private static void transpose(final double[][] source,
                                  final double[][] target,
                                  final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                target[j][i] = source[i][j];
            }
        }
    }

    private static void multiply(final double[][] left,
                                 final double[][] right,
                                 final double[][] out,
                                 final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += left[i][k] * right[k][j];
                }
                out[i][j] = sum;
            }
        }
    }

    private static void multiply(final double[][] left,
                                 final double[][] right,
                                 final double[] out,
                                 final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            final int iOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += left[i][k] * right[k][j];
                }
                out[iOffset + j] = sum;
            }
        }
    }

    private static void multiply(final double[][] matrix,
                                 final double[] vector,
                                 final double[] out,
                                 final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            double sum = 0.0;
            for (int j = 0; j < dimension; ++j) {
                sum += matrix[i][j] * vector[j];
            }
            out[i] = sum;
        }
    }

    private static void copyMatrix(final double[][] source,
                                   final double[] target,
                                   final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            System.arraycopy(source[i], 0, target, i * dimension, dimension);
        }
    }

    private static void scaleInPlace(final double[] vector,
                                     final double factor,
                                     final int count) {
        for (int i = 0; i < count; ++i) {
            vector[i] *= factor;
        }
    }

    private static double dot(final double[] left, final double[] right, final int dimension) {
        double sum = 0.0;
        for (int i = 0; i < dimension; ++i) {
            sum += left[i] * right[i];
        }
        return sum;
    }
}
