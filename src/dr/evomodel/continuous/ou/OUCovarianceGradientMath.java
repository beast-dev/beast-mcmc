package dr.evomodel.continuous.ou;

import java.util.Arrays;

public final class OUCovarianceGradientMath {

    public static final double[] GL5_NODES = {
            0.0,
            -0.5384693101056831, 0.5384693101056831,
            -0.9061798459386640, 0.9061798459386640
    };

    public static final double[] GL5_WEIGHTS = {
            0.5688888888888889,
            0.4786286704993665, 0.4786286704993665,
            0.2369268850561891, 0.2369268850561891
    };

    private OUCovarianceGradientMath() { }

    static void fillSymmetricFromFlat(final double[] source,
                                      final boolean transpose,
                                      final int dim,
                                      final double[] out,
                                      final int outOffset) {
        for (int i = 0; i < dim; ++i) {
            final int rowOffset = outOffset + i * dim;
            for (int j = 0; j < dim; ++j) {
                out[rowOffset + j] = 0.5 * (
                        flatSquareValue(source, i, j, dim, transpose)
                                + flatSquareValue(source, j, i, dim, transpose));
            }
        }
    }

    static double flatSquareValue(final double[] source,
                                  final int row,
                                  final int col,
                                  final int dim,
                                  final boolean transpose) {
        return transpose ? source[col * dim + row] : source[row * dim + col];
    }

    static void clearMatrix(final double[] matrix,
                            final int offset,
                            final int length) {
        Arrays.fill(matrix, offset, offset + length, 0.0);
    }

    /** Computes out = expm(-A * t). */
    public static void buildExpmMinusAs(final double t,
                                        final double[] a,
                                        final int aOffset,
                                        final int dim,
                                        final double[] out,
                                        final int outOffset,
                                        final double[] scaledMinusA,
                                        final int scaledMinusAOffset) {
        for (int i = 0; i < dim; ++i) {
            final int rowOffset = i * dim;
            for (int j = 0; j < dim; ++j) {
                scaledMinusA[scaledMinusAOffset + rowOffset + j] =
                        -t * a[aOffset + rowOffset + j];
            }
        }
        MatrixExponentialUtils.expmFlat(scaledMinusA, scaledMinusAOffset, out, outOffset, dim);
    }

    static void buildVanLoanCovariance(final double t,
                                       final double[] a,
                                       final int aOffset,
                                       final double[] q,
                                       final int qOffset,
                                       final int dim,
                                       final double[] out,
                                       final int outOffset,
                                       final double[] vanLoan,
                                       final int vanLoanOffset,
                                       final double[] exp,
                                       final int expOffset) {
        final int twoDim = 2 * dim;
        for (int i = 0; i < dim; ++i) {
            final int rowOffset = i * dim;
            final int topRow = vanLoanOffset + i * twoDim;
            final int bottomRow = vanLoanOffset + (i + dim) * twoDim;
            for (int j = 0; j < dim; ++j) {
                vanLoan[topRow + j] = -t * a[aOffset + rowOffset + j];
                vanLoan[topRow + j + dim] = t * q[qOffset + rowOffset + j];
                vanLoan[bottomRow + j] = 0.0;
                vanLoan[bottomRow + j + dim] = t * a[aOffset + j * dim + i];
            }
        }
        MatrixExponentialUtils.expmFlat(vanLoan, vanLoanOffset, exp, expOffset, twoDim);
        for (int i = 0; i < dim; ++i) {
            final int outRow = outOffset + i * dim;
            for (int j = 0; j < dim; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dim; ++k) {
                    sum += exp[expOffset + i * twoDim + k + dim]
                            * exp[expOffset + j * twoDim + k];
                }
                out[outRow + j] = sum;
            }
        }
    }

    static double[] allocateMatrixStack(final int count,
                                        final int matrixLength) {
        return new double[count * matrixLength];
    }

    static double[] allocateMatrix(final int rows, final int cols) {
        return new double[rows * cols];
    }
}
