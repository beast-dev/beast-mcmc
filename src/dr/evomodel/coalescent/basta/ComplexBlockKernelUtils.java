package dr.evomodel.coalescent.basta;

import dr.evomodel.substmodel.EigenDecomposition;

final class ComplexBlockKernelUtils {

    private static final double EIGEN_TOLERANCE = 1.0e-12;
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

    private ComplexBlockKernelUtils() { }

    static ComplexKernelPlan buildPlan(EigenDecomposition decomposition, double time, int stateCount) {
        int[] blockStarts = new int[stateCount];
        int[] blockDims = new int[stateCount];
        int blockCount = buildEigenBlocks(decomposition.getEigenValues(), blockStarts, blockDims, stateCount);
        ComplexKernelEntry[] entries = new ComplexKernelEntry[blockCount * blockCount];
        int entryIndex = 0;

        double[] eigenValues = decomposition.getEigenValues();

        for (int leftBlock = 0; leftBlock < blockCount; ++leftBlock) {
            final int leftStart = blockStarts[leftBlock];
            final int leftDim = blockDims[leftBlock];
            for (int rightBlock = 0; rightBlock < blockCount; ++rightBlock) {
                final int rightStart = blockStarts[rightBlock];
                final int rightDim = blockDims[rightBlock];
                double[] coefficients = buildCoefficients(
                        eigenValues, leftStart, leftDim, rightStart, rightDim, stateCount, time);
                entries[entryIndex++] = new ComplexKernelEntry(leftStart, leftDim, rightStart, rightDim, coefficients);
            }
        }
        return new ComplexKernelPlan(entries);
    }

    static void applyPlan(ComplexKernelPlan plan,
                          double[] transformed,
                          double[] eigenBasisGradient,
                          Workspace workspace,
                          int stateCount) {
        for (ComplexKernelEntry entry : plan.entries) {
            int size = entry.leftDim * entry.rightDim;
            double[] in = workspace.kernelInput[size];

            for (int i = 0; i < entry.leftDim; ++i) {
                int transformedRowOffset = (entry.leftStart + i) * stateCount;
                for (int j = 0; j < entry.rightDim; ++j) {
                    in[i * entry.rightDim + j] = transformed[transformedRowOffset + entry.rightStart + j];
                }
            }

            if (size == 1) {
                int gradientOffset = entry.leftStart * stateCount + entry.rightStart;
                eigenBasisGradient[gradientOffset] += entry.coefficients[0] * in[0];
                continue;
            }

            if (size == 2) {
                final double in0 = in[0];
                final double in1 = in[1];
                final double[] c = entry.coefficients;
                final double out0 = c[0] * in0 + c[1] * in1;
                final double out1 = c[2] * in0 + c[3] * in1;
                if (entry.rightDim == 2) {
                    int base = entry.leftStart * stateCount + entry.rightStart;
                    eigenBasisGradient[base] += out0;
                    eigenBasisGradient[base + 1] += out1;
                } else {
                    int base = entry.leftStart * stateCount + entry.rightStart;
                    eigenBasisGradient[base] += out0;
                    eigenBasisGradient[base + stateCount] += out1;
                }
                continue;
            }

            final double in0 = in[0];
            final double in1 = in[1];
            final double in2 = in[2];
            final double in3 = in[3];
            final double[] c = entry.coefficients;
            final double out0 = c[0] * in0 + c[1] * in1 + c[2] * in2 + c[3] * in3;
            final double out1 = c[4] * in0 + c[5] * in1 + c[6] * in2 + c[7] * in3;
            final double out2 = c[8] * in0 + c[9] * in1 + c[10] * in2 + c[11] * in3;
            final double out3 = c[12] * in0 + c[13] * in1 + c[14] * in2 + c[15] * in3;
            int base = entry.leftStart * stateCount + entry.rightStart;
            eigenBasisGradient[base] += out0;
            eigenBasisGradient[base + 1] += out1;
            eigenBasisGradient[base + stateCount] += out2;
            eigenBasisGradient[base + stateCount + 1] += out3;
        }
    }

    static final class ComplexKernelPlan {
        final ComplexKernelEntry[] entries;

        private ComplexKernelPlan(ComplexKernelEntry[] entries) {
            this.entries = entries;
        }
    }

    static final class ComplexKernelEntry {
        final int leftStart;
        final int leftDim;
        final int rightStart;
        final int rightDim;
        final double[] coefficients;

        private ComplexKernelEntry(int leftStart, int leftDim, int rightStart, int rightDim, double[] coefficients) {
            this.leftStart = leftStart;
            this.leftDim = leftDim;
            this.rightStart = rightStart;
            this.rightDim = rightDim;
            this.coefficients = coefficients;
        }
    }

    static final class Workspace {
        final double[][] kernelInput;

        Workspace() {
            this.kernelInput = new double[5][];
            this.kernelInput[1] = new double[1];
            this.kernelInput[2] = new double[2];
            this.kernelInput[4] = new double[4];
        }
    }

    private static int buildEigenBlocks(double[] eigenValues, int[] blockStarts, int[] blockDims, int stateCount) {
        int blockCount = 0;
        for (int i = 0; i < stateCount; ++i) {
            blockStarts[blockCount] = i;
            if (eigenValues.length > stateCount && Math.abs(eigenValues[stateCount + i]) > EIGEN_TOLERANCE) {
                blockDims[blockCount] = 2;
                i++;
            } else {
                blockDims[blockCount] = 1;
            }
            blockCount++;
        }
        return blockCount;
    }

    private static void fillEigenBlock(double[][] localQ,
                                       int offset,
                                       int eigenStart,
                                       int blockDim,
                                       double[] eigenValues,
                                       int stateCount) {
        localQ[offset][offset] = eigenValues[eigenStart];
        if (blockDim == 2) {
            double imag = eigenValues[eigenStart + stateCount];
            localQ[offset][offset + 1] = imag;
            localQ[offset + 1][offset] = -imag;
            localQ[offset + 1][offset + 1] = eigenValues[eigenStart];
        }
    }

    private static void zeroSquare(double[][] matrix, int dim) {
        for (int i = 0; i < dim; ++i) {
            java.util.Arrays.fill(matrix[i], 0, dim, 0.0);
        }
    }

    private static double[] buildCoefficients(double[] eigenValues,
                                              int leftStart,
                                              int leftDim,
                                              int rightStart,
                                              int rightDim,
                                              int stateCount,
                                              double time) {
        if (leftDim == 1 && rightDim == 1) {
            return buildScalarCoefficient(eigenValues[leftStart], eigenValues[rightStart], time);
        }
        if (leftDim == 1 && rightDim == 2) {
            return buildOneByTwoCoefficient(
                    eigenValues[leftStart],
                    eigenValues[rightStart],
                    eigenValues[rightStart + stateCount],
                    time);
        }
        if (leftDim == 2 && rightDim == 1) {
            return buildTwoByOneCoefficient(
                    eigenValues[leftStart],
                    eigenValues[leftStart + stateCount],
                    eigenValues[rightStart],
                    time);
        }
        if (leftDim == 2 && rightDim == 2) {
            return buildTwoByTwoCoefficient(
                    eigenValues[leftStart],
                    eigenValues[leftStart + stateCount],
                    eigenValues[rightStart],
                    eigenValues[rightStart + stateCount],
                    time);
        }

        int size = leftDim * rightDim;
        double[][] aTranspose = new double[leftDim][leftDim];
        double[][] leftExp = new double[leftDim][leftDim];
        double[][] b = new double[rightDim][rightDim];
        fillTransposeEigenBlock(aTranspose, leftExp, leftStart, leftDim, eigenValues, stateCount, time);
        fillOriginalEigenBlock(b, rightStart, rightDim, eigenValues, stateCount);

        double[][] generator = new double[size][size];
        fillRowMajorKroneckerGenerator(generator, aTranspose, b, leftDim, rightDim);

        double[][] phi = phiMatrix(generator, time);
        double[] coefficients = new double[size * size];
        leftMultiplyByExp(coefficients, leftExp, phi, leftDim, rightDim);
        return coefficients;
    }

    private static double[] buildScalarCoefficient(double leftEigenvalue,
                                                   double rightEigenvalue,
                                                   double time) {
        double[] coefficients = new double[1];
        double delta = leftEigenvalue - rightEigenvalue;
        if (Math.abs(delta) < EIGEN_TOLERANCE) {
            coefficients[0] = time * Math.exp(time * leftEigenvalue);
        } else {
            coefficients[0] = (Math.exp(time * leftEigenvalue) - Math.exp(time * rightEigenvalue)) / delta;
        }
        return coefficients;
    }

    private static double[] buildOneByTwoCoefficient(double leftEigenvalue,
                                                     double rightReal,
                                                     double rightImag,
                                                     double time) {
        double[] integral = new double[4];
        fillScaledRotationIntegral(rightReal - leftEigenvalue, rightImag, time, integral);
        double scale = Math.exp(time * leftEigenvalue);
        for (int i = 0; i < integral.length; ++i) {
            integral[i] *= scale;
        }
        return integral;
    }

    private static double[] buildTwoByOneCoefficient(double leftReal,
                                                     double leftImag,
                                                     double rightEigenvalue,
                                                     double time) {
        double[] integral = new double[4];
        fillScaledRotationIntegral(rightEigenvalue - leftReal, leftImag, time, integral);

        double exp = Math.exp(time * leftReal);
        double cos = Math.cos(time * leftImag);
        double sin = Math.sin(time * leftImag);
        double l00 = exp * cos;
        double l01 = -exp * sin;
        double l10 = exp * sin;
        double l11 = exp * cos;

        double[] coefficients = new double[4];
        coefficients[0] = l00 * integral[0] + l01 * integral[2];
        coefficients[1] = l00 * integral[1] + l01 * integral[3];
        coefficients[2] = l10 * integral[0] + l11 * integral[2];
        coefficients[3] = l10 * integral[1] + l11 * integral[3];
        return coefficients;
    }

    private static double[] buildTwoByTwoCoefficient(double leftReal,
                                                     double leftImag,
                                                     double rightReal,
                                                     double rightImag,
                                                     double time) {
        double[] plus = complexMatrix(
                multiplyComplex(
                        expComplex(leftReal, -leftImag, time),
                        integralComplex(rightReal - leftReal, leftImag + rightImag, time)));
        double[] minus = complexMatrix(
                multiplyComplex(
                        expComplex(leftReal, leftImag, time),
                        integralComplex(rightReal - leftReal, rightImag - leftImag, time)));

        double[] coefficients = new double[16];
        double[][] basisColumns = {
                {0.5, 0.0, 0.5, 0.0},
                {0.0, 0.5, 0.0, 0.5},
                {0.0, -0.5, 0.0, 0.5},
                {0.5, 0.0, -0.5, 0.0}
        };

        for (int col = 0; col < 4; ++col) {
            double u = basisColumns[col][0];
            double v = basisColumns[col][1];
            double p = basisColumns[col][2];
            double q = basisColumns[col][3];

            double outU = minus[0] * u + minus[1] * v;
            double outV = minus[2] * u + minus[3] * v;
            double outP = plus[0] * p + plus[1] * q;
            double outQ = plus[2] * p + plus[3] * q;

            coefficients[col] = outU + outP;
            coefficients[4 + col] = outV + outQ;
            coefficients[8 + col] = -outV + outQ;
            coefficients[12 + col] = outU - outP;
        }
        return coefficients;
    }

    private static double[] expComplex(double real, double imag, double time) {
        double scale = Math.exp(time * real);
        return new double[] {scale * Math.cos(time * imag), scale * Math.sin(time * imag)};
    }

    private static double[] integralComplex(double real, double imag, double time) {
        double denom = real * real + imag * imag;
        if (denom < EIGEN_TOLERANCE) {
            return new double[] {time, 0.0};
        }
        double[] exp = expComplex(real, imag, time);
        return new double[] {
                (real * (exp[0] - 1.0) + imag * exp[1]) / denom,
                (real * exp[1] - imag * (exp[0] - 1.0)) / denom
        };
    }

    private static double[] multiplyComplex(double[] left, double[] right) {
        return new double[] {
                left[0] * right[0] - left[1] * right[1],
                left[0] * right[1] + left[1] * right[0]
        };
    }

    private static double[] complexMatrix(double[] z) {
        return new double[] {
                z[0], z[1],
                -z[1], z[0]
        };
    }

    private static void fillScaledRotationIntegral(double realShift,
                                                   double imagShift,
                                                   double time,
                                                   double[] out) {
        double denominator = realShift * realShift + imagShift * imagShift;
        if (denominator < EIGEN_TOLERANCE) {
            out[0] = time;
            out[1] = 0.0;
            out[2] = 0.0;
            out[3] = time;
            return;
        }

        double exp = Math.exp(time * realShift);
        double cos = Math.cos(time * imagShift);
        double sin = Math.sin(time * imagShift);

        double integralCos = (exp * (realShift * cos + imagShift * sin) - realShift) / denominator;
        double integralSin = (exp * (realShift * sin - imagShift * cos) + imagShift) / denominator;

        out[0] = integralCos;
        out[1] = integralSin;
        out[2] = -integralSin;
        out[3] = integralCos;
    }

    private static void fillTransposeEigenBlock(double[][] transposeBlock,
                                                double[][] expTransposeBlock,
                                                int eigenStart,
                                                int blockDim,
                                                double[] eigenValues,
                                                int stateCount,
                                                double time) {
        double real = eigenValues[eigenStart];
        transposeBlock[0][0] = real;
        if (blockDim == 1) {
            expTransposeBlock[0][0] = Math.exp(time * real);
            return;
        }

        double imag = eigenValues[eigenStart + stateCount];
        transposeBlock[0][1] = -imag;
        transposeBlock[1][0] = imag;
        transposeBlock[1][1] = real;

        double exp = Math.exp(time * real);
        double cos = Math.cos(time * imag);
        double sin = Math.sin(time * imag);
        expTransposeBlock[0][0] = exp * cos;
        expTransposeBlock[0][1] = -exp * sin;
        expTransposeBlock[1][0] = exp * sin;
        expTransposeBlock[1][1] = exp * cos;
    }

    private static void fillOriginalEigenBlock(double[][] block,
                                               int eigenStart,
                                               int blockDim,
                                               double[] eigenValues,
                                               int stateCount) {
        double real = eigenValues[eigenStart];
        block[0][0] = real;
        if (blockDim == 1) {
            return;
        }

        double imag = eigenValues[eigenStart + stateCount];
        block[0][1] = imag;
        block[1][0] = -imag;
        block[1][1] = real;
    }

    private static void fillRowMajorKroneckerGenerator(double[][] generator,
                                                       double[][] aTranspose,
                                                       double[][] b,
                                                       int leftDim,
                                                       int rightDim) {
        zeroSquare(generator, generator.length);

        for (int outLeft = 0; outLeft < leftDim; ++outLeft) {
            for (int inLeft = 0; inLeft < leftDim; ++inLeft) {
                double value = -aTranspose[outLeft][inLeft];
                if (value == 0.0) {
                    continue;
                }
                int outBase = outLeft * rightDim;
                int inBase = inLeft * rightDim;
                for (int right = 0; right < rightDim; ++right) {
                    generator[outBase + right][inBase + right] += value;
                }
            }
        }

        for (int left = 0; left < leftDim; ++left) {
            int base = left * rightDim;
            for (int outRight = 0; outRight < rightDim; ++outRight) {
                int outIndex = base + outRight;
                for (int inRight = 0; inRight < rightDim; ++inRight) {
                    generator[outIndex][base + inRight] += b[outRight][inRight];
                }
            }
        }
    }

    private static double[][] phiMatrix(double[][] generator, double time) {
        int size = generator.length;
        int blockSize = 2 * size;
        double[][] block = new double[blockSize][blockSize];
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                block[i][j] = time * generator[i][j];
            }
            block[i][size + i] = 1.0;
        }

        double[][] expBlock = expmPade13(block);
        double[][] phi = new double[size][size];
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                phi[i][j] = time * expBlock[i][size + j];
            }
        }
        return phi;
    }

    private static void leftMultiplyByExp(double[] coefficients,
                                          double[][] leftExp,
                                          double[][] phi,
                                          int leftDim,
                                          int rightDim) {
        int size = leftDim * rightDim;
        if (leftDim == 1) {
            double scale = leftExp[0][0];
            for (int i = 0; i < size; ++i) {
                int offset = i * size;
                for (int j = 0; j < size; ++j) {
                    coefficients[offset + j] = scale * phi[i][j];
                }
            }
            return;
        }

        for (int outLeft = 0; outLeft < leftDim; ++outLeft) {
            int outBase = outLeft * rightDim;
            for (int inLeft = 0; inLeft < leftDim; ++inLeft) {
                double scale = leftExp[outLeft][inLeft];
                if (scale == 0.0) {
                    continue;
                }
                int midBase = inLeft * rightDim;
                for (int right = 0; right < rightDim; ++right) {
                    int outRow = outBase + right;
                    int phiRow = midBase + right;
                    int coeffOffset = outRow * size;
                    for (int col = 0; col < size; ++col) {
                        coefficients[coeffOffset + col] += scale * phi[phiRow][col];
                    }
                }
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
        for (int i = 0; i < scale; ++i) {
            result = multiply(result, result);
        }
        return result;
    }

    private static double matrixOneNorm(double[][] matrix) {
        double max = 0.0;
        for (int j = 0; j < matrix[0].length; j++) {
            double sum = 0.0;
            for (double[] row : matrix) {
                sum += Math.abs(row[j]);
            }
            max = Math.max(max, sum);
        }
        return max;
    }

    private static double[][] identity(int n) {
        double[][] out = new double[n][n];
        for (int i = 0; i < n; ++i) {
            out[i][i] = 1.0;
        }
        return out;
    }

    private static double[][] scale(double[][] matrix, double factor) {
        int n = matrix.length;
        int m = matrix[0].length;
        double[][] out = new double[n][m];
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
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
        for (int i = 0; i < n; ++i) {
            for (int k = 0; k < inner; ++k) {
                double aik = a[i][k];
                if (aik == 0.0) {
                    continue;
                }
                for (int j = 0; j < m; ++j) {
                    out[i][j] += aik * b[k][j];
                }
            }
        }
        return out;
    }

    private static double[][] add(double[][] a, double[][] b) {
        double[][] out = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                out[i][j] = a[i][j] + b[i][j];
            }
        }
        return out;
    }

    private static double[][] subtract(double[][] a, double[][] b) {
        double[][] out = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                out[i][j] = a[i][j] - b[i][j];
            }
        }
        return out;
    }

    private static double[][] solve(double[][] a, double[][] b) {
        int n = a.length;
        int m = b[0].length;
        double[][] aug = new double[n][n + m];
        for (int i = 0; i < n; ++i) {
            System.arraycopy(a[i], 0, aug[i], 0, n);
            System.arraycopy(b[i], 0, aug[i], n, m);
        }

        for (int col = 0; col < n; ++col) {
            int pivot = col;
            double max = Math.abs(aug[col][col]);
            for (int row = col + 1; row < n; ++row) {
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
            for (int j = col; j < n + m; ++j) {
                aug[col][j] /= pivotValue;
            }

            for (int row = 0; row < n; ++row) {
                if (row == col) {
                    continue;
                }
                double factor = aug[row][col];
                if (factor == 0.0) {
                    continue;
                }
                for (int j = col; j < n + m; ++j) {
                    aug[row][j] -= factor * aug[col][j];
                }
            }
        }

        double[][] out = new double[n][m];
        for (int i = 0; i < n; ++i) {
            System.arraycopy(aug[i], n, out[i], 0, m);
        }
        return out;
    }

    private static double log2(double x) {
        return Math.log(x) / Math.log(2.0);
    }
}
