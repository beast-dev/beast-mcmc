package dr.evomodel.coalescent.basta;

final class AdjointMatrixExponentialUtils {

    private static final double[] PADE13 = new double[] {
            64764752532480000.0,
            32382376266240000.0,
            7771770303897600.0,
            1187353796428800.0,
            129060195264000.0,
            10559470521600.0,
            670442572800.0,
            33522128640.0,
            1323241920.0,
            40840800.0,
            960960.0,
            16380.0,
            182.0,
            1.0
    };

    private static final double THETA13 = 5.371920351148152;

    private AdjointMatrixExponentialUtils() { }

    static void accumulateExpmAdjoint(double[][] q,
                                      double time,
                                      double[][] matrixAdjoint,
                                      double[][] rateGradient,
                                      Workspace workspace) {
        buildFrechetBlockFromTransposeScaled(q, time, matrixAdjoint, workspace.block);
        double[][] expBlock = expmPade13(workspace.block);
        int n = q.length;
        for (int i = 0; i < n; i++) {
            System.arraycopy(expBlock[i], n, workspace.frechet[i], 0, n);
        }
        addScaledInPlace(rateGradient, workspace.frechet, time);
    }

    static final class Workspace {
        final double[][] block;
        final double[][] frechet;

        Workspace(int stateCount) {
            this.block = new double[2 * stateCount][2 * stateCount];
            this.frechet = new double[stateCount][stateCount];
        }
    }

    private static void buildFrechetBlockFromTransposeScaled(double[][] q,
                                                             double time,
                                                             double[][] direction,
                                                             double[][] block) {
        int n = q.length;
        clear(block);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double scaledTranspose = time * q[j][i];
                block[i][j] = scaledTranspose;
                block[i][n + j] = direction[i][j];
                block[n + i][n + j] = scaledTranspose;
            }
        }
    }

    private static double[][] expmPade13(double[][] matrix) {
        int n = matrix.length;
        double[][] ident = identity(n);
        double norm = matrixOneNorm(matrix);
        int scale = norm > 0.0 ? Math.max(0, (int) Math.ceil(log2(norm / THETA13))) : 0;
        double[][] a = scale(matrix, 1.0 / Math.pow(2.0, scale));

        double[][] a2 = multiply(a, a);
        double[][] a4 = multiply(a2, a2);
        double[][] a6 = multiply(a4, a2);

        double[][] uInner = add(
                add(multiply(a6, add(add(scale(a6, PADE13[13]), scale(a4, PADE13[11])), scale(a2, PADE13[9]))), scale(a6, PADE13[7])),
                add(scale(a4, PADE13[5]), add(scale(a2, PADE13[3]), scale(ident, PADE13[1])))
        );
        double[][] u = multiply(a, uInner);

        double[][] v = add(
                add(multiply(a6, add(add(scale(a6, PADE13[12]), scale(a4, PADE13[10])), scale(a2, PADE13[8]))), scale(a6, PADE13[6])),
                add(scale(a4, PADE13[4]), add(scale(a2, PADE13[2]), scale(ident, PADE13[0])))
        );

        double[][] result = solve(subtract(v, u), add(v, u));
        for (int i = 0; i < scale; i++) {
            result = multiply(result, result);
        }
        return result;
    }

    private static void clear(double[][] matrix) {
        for (double[] row : matrix) {
            java.util.Arrays.fill(row, 0.0);
        }
    }

    private static double matrixOneNorm(double[][] matrix) {
        double max = 0.0;
        for (int j = 0; j < matrix[0].length; j++) {
            double sum = 0.0;
            for (double[] doubles : matrix) {
                sum += Math.abs(doubles[j]);
            }
            max = Math.max(max, sum);
        }
        return max;
    }

    private static double[][] identity(int n) {
        double[][] out = new double[n][n];
        for (int i = 0; i < n; i++) {
            out[i][i] = 1.0;
        }
        return out;
    }

    private static double[][] scale(double[][] matrix, double factor) {
        int n = matrix.length;
        int m = matrix[0].length;
        double[][] out = new double[n][m];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                out[i][j] = factor * matrix[i][j];
            }
        }
        return out;
    }

    private static double[][] multiply(double[][] a, double[][] b) {
        int n = a.length;
        int m = b[0].length;
        int inner = b.length;
        double[][] out = new double[n][m];
        for (int i = 0; i < n; i++) {
            for (int k = 0; k < inner; k++) {
                double aik = a[i][k];
                if (aik == 0.0) {
                    continue;
                }
                for (int j = 0; j < m; j++) {
                    out[i][j] += aik * b[k][j];
                }
            }
        }
        return out;
    }

    private static double[][] add(double[][] a, double[][] b) {
        double[][] out = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                out[i][j] = a[i][j] + b[i][j];
            }
        }
        return out;
    }

    private static double[][] subtract(double[][] a, double[][] b) {
        double[][] out = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                out[i][j] = a[i][j] - b[i][j];
            }
        }
        return out;
    }

    private static void addScaledInPlace(double[][] target, double[][] source, double factor) {
        for (int i = 0; i < target.length; i++) {
            for (int j = 0; j < target[i].length; j++) {
                target[i][j] += factor * source[i][j];
            }
        }
    }

    private static double[][] solve(double[][] a, double[][] b) {
        int n = a.length;
        int m = b[0].length;
        double[][] aug = new double[n][n + m];
        for (int i = 0; i < n; i++) {
            System.arraycopy(a[i], 0, aug[i], 0, n);
            System.arraycopy(b[i], 0, aug[i], n, m);
        }

        for (int col = 0; col < n; col++) {
            int pivot = col;
            double max = Math.abs(aug[col][col]);
            for (int row = col + 1; row < n; row++) {
                double value = Math.abs(aug[row][col]);
                if (value > max) {
                    max = value;
                    pivot = row;
                }
            }
            if (pivot != col) {
                double[] tmp = aug[col];
                aug[col] = aug[pivot];
                aug[pivot] = tmp;
            }

            double pivotValue = aug[col][col];
            for (int j = col; j < n + m; j++) {
                aug[col][j] /= pivotValue;
            }

            for (int row = 0; row < n; row++) {
                if (row == col) {
                    continue;
                }
                double factor = aug[row][col];
                if (factor == 0.0) {
                    continue;
                }
                for (int j = col; j < n + m; j++) {
                    aug[row][j] -= factor * aug[col][j];
                }
            }
        }

        double[][] out = new double[n][m];
        for (int i = 0; i < n; i++) {
            System.arraycopy(aug[i], n, out[i], 0, m);
        }
        return out;
    }

    private static double log2(double x) {
        return Math.log(x) / Math.log(2.0);
    }
}
