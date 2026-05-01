package dr.evomodel.continuous.ou;

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

    static void fillSymmetricFromMatrix(final double[][] source,
                                        final int dim,
                                        final double[][] out) {
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                out[i][j] = 0.5 * (source[i][j] + source[j][i]);
            }
        }
    }

    static void fillSymmetricFromFlat(final double[] source,
                                      final boolean transpose,
                                      final int dim,
                                      final double[][] out) {
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                out[i][j] = 0.5 * (
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

    static void clearSquare(final double[][] matrix) {
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = 0; j < matrix[i].length; ++j) {
                matrix[i][j] = 0.0;
            }
        }
    }

    /** Computes out = expm(-A * t). */
    public static void buildExpmMinusAs(final double t,
                                        final double[][] a,
                                        final int dim,
                                        final double[][] out,
                                        final double[][] scaledMinusA) {
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                scaledMinusA[i][j] = -t * a[i][j];
            }
        }
        MatrixExponentialUtils.expm(scaledMinusA, out);
    }

    static void buildVanLoanCovariance(final double t,
                                       final double[][] a,
                                       final double[][] q,
                                       final int dim,
                                       final double[][] out,
                                       final double[][] vanLoan,
                                       final double[][] exp) {
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                vanLoan[i][j] = -t * a[i][j];
                vanLoan[i][j + dim] = t * q[i][j];
                vanLoan[i + dim][j] = 0.0;
                vanLoan[i + dim][j + dim] = t * a[j][i];
            }
        }
        MatrixExponentialUtils.expm(vanLoan, exp);
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dim; ++k) {
                    sum += exp[i][k + dim] * exp[j][k];
                }
                out[i][j] = sum;
            }
        }
    }

    static double[][][] allocateMatrixStack(final int count,
                                            final int rows,
                                            final int cols) {
        final double[][][] stack = new double[count][][];
        for (int i = 0; i < count; ++i) {
            stack[i] = allocateMatrix(rows, cols);
        }
        return stack;
    }

    static double[][] allocateMatrix(final int rows, final int cols) {
        final double[][] matrix = new double[rows][];
        for (int i = 0; i < rows; ++i) {
            matrix[i] = new double[cols];
        }
        return matrix;
    }
}
