package dr.evomodel.treedatalikelihood.continuous.backprop;

import org.ejml.data.DenseMatrix64F;

/**
 * Exact cached coefficient plan for the block Fréchet map
 *
 * <pre>
 *   E -> ∫₀¹ exp((1-s)A) E exp(sB) ds
 * </pre>
 *
 * on all 1x1 / 2x2 block pairs of a block-diagonal matrix.
 *
 * <p>The production path here is fully closed form:
 * 1x1 / 1x2 / 2x1 pairs are evaluated through analytic block phi functions,
 * and 2x2 / 2x2 pairs are assembled from the exact Sylvester identity using
 * explicit 4x4 linear solves. There is no numerical quadrature and no generic
 * matrix exponential in this class.</p>
 */
final class BlockDiagonalFrechetExactPlan {

    private static final double EPS = 1.0e-12;
    private static final double SMALL_ROOT_SERIES_CUTOFF = 1.0e-3;
    private static final int SMALL_ROOT_SERIES_ORDER = 4;
    private static final double MOMENT_SERIES_CUTOFF = 1.0e-7;
    private static final int MOMENT_SERIES_TERMS = 30;
    private static final int MAX_SMALL_ROOT_DERIVATIVE_ORDER = 2 * SMALL_ROOT_SERIES_ORDER + 1;
    private static final int MAX_SMALL_ROOT_MOMENT_POWER = 2 * MAX_SMALL_ROOT_DERIVATIVE_ORDER;
    private static final int MAX_FACTORIAL_INDEX = 2 * SMALL_ROOT_SERIES_ORDER + 1;
    private static final double[] INVERSE_FACTORIALS = buildInverseFactorials(MAX_FACTORIAL_INDEX);
    private static final double[][] BINOMIAL = buildBinomialTable(MAX_SMALL_ROOT_DERIVATIVE_ORDER);
    private static final ThreadLocal<ComplexWorkspace> COMPLEX_WORKSPACE =
            ThreadLocal.withInitial(ComplexWorkspace::new);
    private static final ThreadLocal<RealSeriesWorkspace> REAL_SERIES_WORKSPACE =
            ThreadLocal.withInitial(RealSeriesWorkspace::new);

    private BlockDiagonalFrechetExactPlan() {
    }

    static Plan buildPlan(final BlockDiagonalExpSolver.BlockStructure structure,
                          final double[] blockDParams,
                          final double t) {
        final int dim = structure.getDim();
        final int blockCount = structure.getNumBlocks();
        final Entry[] entries = new Entry[blockCount * blockCount];
        int entryIndex = 0;

        for (int leftBlock = 0; leftBlock < blockCount; ++leftBlock) {
            final int leftStart = structure.getBlockStart(leftBlock);
            final int leftDim = structure.getBlockSize(leftBlock);
            for (int rightBlock = 0; rightBlock < blockCount; ++rightBlock) {
                final int rightStart = structure.getBlockStart(rightBlock);
                final int rightDim = structure.getBlockSize(rightBlock);
                entries[entryIndex++] = new Entry(
                        leftStart,
                        leftDim,
                        rightStart,
                        rightDim,
                        buildCoefficients(blockDParams, dim, leftStart, leftDim, rightStart, rightDim, t));
            }
        }
        return new Plan(entries);
    }

    static final class Plan {
        private final Entry[] entries;

        private Plan(final Entry[] entries) {
            this.entries = entries;
        }

        void apply(final DenseMatrix64F input,
                   final DenseMatrix64F out,
                   final Workspace workspace) {
            final int dimension = input.numCols;
            final double[] inData = input.data;
            final double[] outData = out.data;

            for (final Entry entry : entries) {
                final int size = entry.size;
                final double[] in = workspace.kernelInput[size];
                final double[] outBlock = workspace.kernelOutput[size];

                for (int i = 0; i < entry.leftDim; ++i) {
                    final int inputBase = (entry.leftStart + i) * dimension + entry.rightStart;
                    for (int j = 0; j < entry.rightDim; ++j) {
                        in[i * entry.rightDim + j] = inData[inputBase + j];
                    }
                }

                if (size == 1) {
                    outBlock[0] = entry.coefficients[0] * in[0];
                } else {
                    for (int row = 0; row < size; ++row) {
                        double sum = 0.0;
                        final int coeffOffset = row * size;
                        for (int col = 0; col < size; ++col) {
                            sum += entry.coefficients[coeffOffset + col] * in[col];
                        }
                        outBlock[row] = sum;
                    }
                }

                for (int i = 0; i < entry.leftDim; ++i) {
                    final int outputBase = (entry.leftStart + i) * dimension + entry.rightStart;
                    for (int j = 0; j < entry.rightDim; ++j) {
                        outData[outputBase + j] = outBlock[i * entry.rightDim + j];
                    }
                }
            }
        }
    }

    static final class Workspace {
        final double[][] kernelInput;
        final double[][] kernelOutput;

        Workspace() {
            this.kernelInput = new double[5][];
            this.kernelOutput = new double[5][];
            this.kernelInput[1] = new double[1];
            this.kernelInput[2] = new double[2];
            this.kernelInput[4] = new double[4];
            this.kernelOutput[1] = new double[1];
            this.kernelOutput[2] = new double[2];
            this.kernelOutput[4] = new double[4];
        }
    }

    private static final class Entry {
        private final int leftStart;
        private final int leftDim;
        private final int rightStart;
        private final int rightDim;
        private final int size;
        private final double[] coefficients;

        private Entry(final int leftStart,
                      final int leftDim,
                      final int rightStart,
                      final int rightDim,
                      final double[] coefficients) {
            this.leftStart = leftStart;
            this.leftDim = leftDim;
            this.rightStart = rightStart;
            this.rightDim = rightDim;
            this.size = leftDim * rightDim;
            this.coefficients = coefficients;
        }
    }

    private static double[] buildCoefficients(final double[] blockDParams,
                                              final int dimension,
                                              final int leftStart,
                                              final int leftDim,
                                              final int rightStart,
                                              final int rightDim,
                                              final double t) {
        final double[] left = new double[4];
        final double[] right = new double[4];
        final double[] leftExp = new double[4];
        final double[] rightExp = new double[4];

        fillScaledBlock(blockDParams, dimension, leftStart, leftDim, -t, left);
        fillScaledBlock(blockDParams, dimension, rightStart, rightDim, -t, right);
        fillSmallExponential(left, leftDim, leftExp);
        fillSmallExponential(right, rightDim, rightExp);

        if (leftDim == 1 && rightDim == 1) {
            return new double[]{dividedDifferenceExp(left[0], right[0])};
        }
        if (leftDim == 1 && rightDim == 2) {
            return buildOneByTwoCoefficient(left[0], right);
        }
        if (leftDim == 2 && rightDim == 1) {
            return buildTwoByOneCoefficient(left, right[0]);
        }
        if (leftDim == 2 && rightDim == 2) {
            return buildTwoByTwoCoefficient(left, right, leftExp, rightExp);
        }

        throw new IllegalStateException("Unsupported block pair: " + leftDim + "x" + rightDim);
    }

    private static void fillScaledBlock(final double[] blockDParams,
                                        final int dimension,
                                        final int start,
                                        final int blockDim,
                                        final double scale,
                                        final double[] out) {
        out[0] = scale * blockDParams[start];
        if (blockDim == 1) {
            return;
        }

        final int upperOffset = dimension;
        final int lowerOffset = dimension + (dimension - 1);
        out[1] = scale * blockDParams[upperOffset + start];
        out[2] = scale * blockDParams[lowerOffset + start];
        out[3] = scale * blockDParams[start + 1];
    }

    private static void fillSmallExponential(final double[] block,
                                             final int blockDim,
                                             final double[] out) {
        if (blockDim == 1) {
            out[0] = Math.exp(block[0]);
            return;
        }

        BlockDiagonalExpSolver.fillExp2x2GeneralBlock(
                -block[0], -block[1], -block[2], -block[3],
                1.0, out, 2, 0);
    }

    private static double[] buildOneByTwoCoefficient(final double leftScalar,
                                                     final double[] rightBlock) {
        final double[] shifted = new double[]{
                rightBlock[0] - leftScalar,
                rightBlock[1],
                rightBlock[2],
                rightBlock[3] - leftScalar
        };
        final double[] phi = new double[4];
        fillPhi2x2(shifted, phi);

        final double scale = Math.exp(leftScalar);
        phi[0] *= scale;
        phi[1] *= scale;
        phi[2] *= scale;
        phi[3] *= scale;

        return new double[]{
                phi[0], phi[2],
                phi[1], phi[3]
        };
    }

    private static double[] buildTwoByOneCoefficient(final double[] leftBlock,
                                                     final double rightScalar) {
        final double[] shifted = new double[]{
                leftBlock[0] - rightScalar,
                leftBlock[1],
                leftBlock[2],
                leftBlock[3] - rightScalar
        };
        final double[] phi = new double[4];
        fillPhi2x2(shifted, phi);

        final double scale = Math.exp(rightScalar);
        for (int i = 0; i < 4; ++i) {
            phi[i] *= scale;
        }
        return phi;
    }

    private static double[] buildTwoByTwoCoefficient(final double[] leftBlock,
                                                     final double[] rightBlock,
                                                     final double[] leftExp,
                                                     final double[] rightExp) {
        if (hasEqualDiagonals(leftBlock) && hasEqualDiagonals(rightBlock)) {
            return buildTwoByTwoCoefficientEqualDiagonal(leftBlock, rightBlock);
        }
        return buildTwoByTwoCoefficientGeneric(leftBlock, rightBlock, leftExp, rightExp);
    }

    private static double[] buildTwoByTwoCoefficientEqualDiagonal(final double[] leftBlock,
                                                                  final double[] rightBlock) {
        final double a = leftBlock[0];
        final double b = rightBlock[0];
        final double leftUpper = leftBlock[1];
        final double leftLower = leftBlock[2];
        final double rightUpper = rightBlock[1];
        final double rightLower = rightBlock[2];

        final double leftProduct = leftUpper * leftLower;
        final double rightProduct = rightUpper * rightLower;
        final double rRe = leftProduct >= 0.0 ? Math.sqrt(leftProduct) : 0.0;
        final double rIm = leftProduct >= 0.0 ? 0.0 : Math.sqrt(-leftProduct);
        final double qRe = rightProduct >= 0.0 ? Math.sqrt(rightProduct) : 0.0;
        final double qIm = rightProduct >= 0.0 ? 0.0 : Math.sqrt(-rightProduct);
        final double rAbs = Math.sqrt(Math.abs(leftProduct));
        final double qAbs = Math.sqrt(Math.abs(rightProduct));
        final ComplexWorkspace workspace = COMPLEX_WORKSPACE.get();

        final EqualDiagonalCoefficients prepared;
        if (rAbs < SMALL_ROOT_SERIES_CUTOFF && qAbs < SMALL_ROOT_SERIES_CUTOFF) {
            prepared = coefficientsBivariateSmallRootSeries(a, b, leftProduct, rightProduct);
        } else if (rAbs < SMALL_ROOT_SERIES_CUTOFF) {
            prepared = coefficientsSmallLeftRootSeries(a, b, qRe, qIm, leftProduct, workspace);
        } else if (qAbs < SMALL_ROOT_SERIES_CUTOFF) {
            prepared = coefficientsSmallRightRootSeries(a, b, rRe, rIm, rightProduct, workspace);
        } else {
            prepared = coefficientsDistinct(a, b, rRe, rIm, qRe, qIm, workspace);
        }

        return buildTwoByTwoCoefficientFromPrepared(prepared, leftUpper, leftLower, rightUpper, rightLower);
    }

    private static double[] buildTwoByTwoCoefficientGeneric(final double[] leftBlock,
                                                            final double[] rightBlock,
                                                            final double[] leftExp,
                                                            final double[] rightExp) {
        try {
            final double[] sylvester = new double[16];
            fillTwoByTwoSylvesterOperator(leftBlock, rightBlock, sylvester);

            final double[] coefficients = new double[16];
            final double[] rhs = new double[4];
            final double[] solution = new double[4];
            final double[] basisMatrix = new double[4];

            for (int basis = 0; basis < 4; ++basis) {
                fillTwoByTwoBasisMatrix(basis, basisMatrix);
                fillTwoByTwoRhs(leftExp, rightExp, basisMatrix, rhs);
                solve4x4(sylvester, rhs, solution);
                coefficients[basis] = solution[0];
                coefficients[4 + basis] = solution[1];
                coefficients[8 + basis] = solution[2];
                coefficients[12 + basis] = solution[3];
            }

            return coefficients;
        } catch (final IllegalStateException degenerate) {
            final Spectrum2x2 leftSpectrum = Spectrum2x2.decompose(leftBlock);
            final Spectrum2x2 rightSpectrum = Spectrum2x2.decompose(rightBlock);
            if (leftSpectrum.distinct && rightSpectrum.distinct) {
                return buildTwoByTwoCoefficientSpectral(leftBlock, rightBlock, leftSpectrum, rightSpectrum);
            }
            if (matricesEqual(leftBlock, rightBlock)) {
                return buildTwoByTwoCoefficientRepeatedSame(leftBlock);
            }
            throw new IllegalStateException(
                    "Degenerate 2x2 Sylvester block pair left=[" + leftBlock[0] + ", " + leftBlock[1] + ", "
                            + leftBlock[2] + ", " + leftBlock[3] + "] right=["
                            + rightBlock[0] + ", " + rightBlock[1] + ", "
                            + rightBlock[2] + ", " + rightBlock[3] + "] leftDistinct="
                            + leftSpectrum.distinct + " rightDistinct=" + rightSpectrum.distinct
                            + " equal=" + matricesEqual(leftBlock, rightBlock),
                    degenerate);
        }
    }

    private static double[] buildTwoByTwoCoefficientFromPrepared(final EqualDiagonalCoefficients prepared,
                                                                 final double leftUpper,
                                                                 final double leftLower,
                                                                 final double rightUpper,
                                                                 final double rightLower) {
        final double[] coefficients = new double[16];
        final double[] basis = new double[4];
        for (int basisIndex = 0; basisIndex < 4; ++basisIndex) {
            fillTwoByTwoBasisMatrix(basisIndex, basis);
            coefficients[basisIndex] =
                    prepared.alpha * basis[0]
                            + prepared.beta * (leftUpper * basis[2])
                            + prepared.gamma * (basis[1] * rightLower)
                            + prepared.eta * (leftUpper * basis[3] * rightLower);
            coefficients[4 + basisIndex] =
                    prepared.alpha * basis[1]
                            + prepared.beta * (leftUpper * basis[3])
                            + prepared.gamma * (basis[0] * rightUpper)
                            + prepared.eta * (leftUpper * basis[2] * rightUpper);
            coefficients[8 + basisIndex] =
                    prepared.alpha * basis[2]
                            + prepared.beta * (leftLower * basis[0])
                            + prepared.gamma * (basis[3] * rightLower)
                            + prepared.eta * (leftLower * basis[1] * rightLower);
            coefficients[12 + basisIndex] =
                    prepared.alpha * basis[3]
                            + prepared.beta * (leftLower * basis[1])
                            + prepared.gamma * (basis[2] * rightUpper)
                            + prepared.eta * (leftLower * basis[0] * rightUpper);
        }
        return coefficients;
    }

    private static EqualDiagonalCoefficients coefficientsDistinct(final double a,
                                                                  final double b,
                                                                  final double rRe,
                                                                  final double rIm,
                                                                  final double qRe,
                                                                  final double qIm,
                                                                  final ComplexWorkspace workspace) {
        final MutableComplex phiPP = workspace.c0;
        final MutableComplex phiPM = workspace.c1;
        final MutableComplex phiMP = workspace.c2;
        final MutableComplex phiMM = workspace.c3;
        final MutableComplex scratch = workspace.c4;

        fillPhiIntegral(a + rRe, rIm, b + qRe, qIm, phiPP, workspace);
        fillPhiIntegral(a + rRe, rIm, b - qRe, -qIm, phiPM, workspace);
        fillPhiIntegral(a - rRe, -rIm, b + qRe, qIm, phiMP, workspace);
        fillPhiIntegral(a - rRe, -rIm, b - qRe, -qIm, phiMM, workspace);

        final double alpha = 0.25 * (phiPP.re + phiPM.re + phiMP.re + phiMM.re);

        divideComplex(
                phiPP.re + phiPM.re - phiMP.re - phiMM.re,
                phiPP.im + phiPM.im - phiMP.im - phiMM.im,
                rRe, rIm, scratch);
        final double beta = 0.25 * scratch.re;

        divideComplex(
                phiPP.re - phiPM.re + phiMP.re - phiMM.re,
                phiPP.im - phiPM.im + phiMP.im - phiMM.im,
                qRe, qIm, scratch);
        final double gamma = 0.25 * scratch.re;

        multiplyComplex(rRe, rIm, qRe, qIm, scratch);
        divideComplex(
                phiPP.re - phiPM.re - phiMP.re + phiMM.re,
                phiPP.im - phiPM.im - phiMP.im + phiMM.im,
                scratch.re, scratch.im, scratch);
        final double eta = 0.25 * scratch.re;

        return new EqualDiagonalCoefficients(alpha, beta, gamma, eta);
    }

    private static EqualDiagonalCoefficients coefficientsBivariateSmallRootSeries(final double a,
                                                                                  final double b,
                                                                                  final double rSquared,
                                                                                  final double qSquared) {
        final RealSeriesWorkspace workspace = REAL_SERIES_WORKSPACE.get();
        fillRealMomentTable(b - a, workspace.moments);
        fillRealDerivativeTable(Math.exp(a), workspace.moments, workspace.derivatives);

        double alpha = 0.0;
        double beta = 0.0;
        double gamma = 0.0;
        double eta = 0.0;

        double rPower = 1.0;
        for (int i = 0; i <= SMALL_ROOT_SERIES_ORDER; ++i) {
            final double rEven = rPower * INVERSE_FACTORIALS[2 * i];
            final double rOdd = rPower * INVERSE_FACTORIALS[2 * i + 1];

            double qPower = 1.0;
            for (int j = 0; j <= SMALL_ROOT_SERIES_ORDER; ++j) {
                final double qEven = qPower * INVERSE_FACTORIALS[2 * j];
                final double qOdd = qPower * INVERSE_FACTORIALS[2 * j + 1];

                alpha += rEven * qEven * workspace.derivatives[derivativeIndex(2 * i, 2 * j)];
                beta += rOdd * qEven * workspace.derivatives[derivativeIndex(2 * i + 1, 2 * j)];
                gamma += rEven * qOdd * workspace.derivatives[derivativeIndex(2 * i, 2 * j + 1)];
                eta += rOdd * qOdd * workspace.derivatives[derivativeIndex(2 * i + 1, 2 * j + 1)];

                qPower *= qSquared;
            }
            rPower *= rSquared;
        }

        return new EqualDiagonalCoefficients(alpha, beta, gamma, eta);
    }

    private static EqualDiagonalCoefficients coefficientsSmallLeftRootSeries(final double a,
                                                                             final double b,
                                                                             final double qRe,
                                                                             final double qIm,
                                                                             final double rSquared,
                                                                             final ComplexWorkspace workspace) {
        double alpha = 0.0;
        double beta = 0.0;
        double gamma = 0.0;
        double eta = 0.0;

        final MutableComplex plus = workspace.c0;
        final MutableComplex minus = workspace.c1;
        final MutableComplex scratch = workspace.c2;

        double rPower = 1.0;
        for (int i = 0; i <= SMALL_ROOT_SERIES_ORDER; ++i) {
            fillMixedDerivative(a, 0.0, b + qRe, qIm, 2 * i, 0, plus, workspace);
            fillMixedDerivative(a, 0.0, b - qRe, -qIm, 2 * i, 0, minus, workspace);
            final double evenAverageRe = 0.5 * (plus.re + minus.re);
            divideComplex(plus.re - minus.re, plus.im - minus.im, qRe, qIm, scratch);
            final double oddQEvenR = 0.5 * scratch.re;

            fillMixedDerivative(a, 0.0, b + qRe, qIm, 2 * i + 1, 0, plus, workspace);
            fillMixedDerivative(a, 0.0, b - qRe, -qIm, 2 * i + 1, 0, minus, workspace);
            final double oddAverageRe = 0.5 * (plus.re + minus.re);
            divideComplex(plus.re - minus.re, plus.im - minus.im, qRe, qIm, scratch);
            final double oddQOddR = 0.5 * scratch.re;

            alpha += rPower * INVERSE_FACTORIALS[2 * i] * evenAverageRe;
            beta += rPower * INVERSE_FACTORIALS[2 * i + 1] * oddAverageRe;
            gamma += rPower * INVERSE_FACTORIALS[2 * i] * oddQEvenR;
            eta += rPower * INVERSE_FACTORIALS[2 * i + 1] * oddQOddR;
            rPower *= rSquared;
        }

        return new EqualDiagonalCoefficients(alpha, beta, gamma, eta);
    }

    private static EqualDiagonalCoefficients coefficientsSmallRightRootSeries(final double a,
                                                                              final double b,
                                                                              final double rRe,
                                                                              final double rIm,
                                                                              final double qSquared,
                                                                              final ComplexWorkspace workspace) {
        double alpha = 0.0;
        double beta = 0.0;
        double gamma = 0.0;
        double eta = 0.0;

        final MutableComplex plus = workspace.c0;
        final MutableComplex minus = workspace.c1;
        final MutableComplex scratch = workspace.c2;

        double qPower = 1.0;
        for (int j = 0; j <= SMALL_ROOT_SERIES_ORDER; ++j) {
            fillMixedDerivative(a + rRe, rIm, b, 0.0, 0, 2 * j, plus, workspace);
            fillMixedDerivative(a - rRe, -rIm, b, 0.0, 0, 2 * j, minus, workspace);
            final double evenAverageRe = 0.5 * (plus.re + minus.re);
            divideComplex(plus.re - minus.re, plus.im - minus.im, rRe, rIm, scratch);
            final double oddREvenQ = 0.5 * scratch.re;

            fillMixedDerivative(a + rRe, rIm, b, 0.0, 0, 2 * j + 1, plus, workspace);
            fillMixedDerivative(a - rRe, -rIm, b, 0.0, 0, 2 * j + 1, minus, workspace);
            final double oddAverageRe = 0.5 * (plus.re + minus.re);
            divideComplex(plus.re - minus.re, plus.im - minus.im, rRe, rIm, scratch);
            final double oddROddQ = 0.5 * scratch.re;

            alpha += qPower * INVERSE_FACTORIALS[2 * j] * evenAverageRe;
            beta += qPower * INVERSE_FACTORIALS[2 * j] * oddREvenQ;
            gamma += qPower * INVERSE_FACTORIALS[2 * j + 1] * oddAverageRe;
            eta += qPower * INVERSE_FACTORIALS[2 * j + 1] * oddROddQ;
            qPower *= qSquared;
        }

        return new EqualDiagonalCoefficients(alpha, beta, gamma, eta);
    }

    private static int derivativeIndex(final int leftOrder,
                                       final int rightOrder) {
        return leftOrder * (MAX_SMALL_ROOT_DERIVATIVE_ORDER + 1) + rightOrder;
    }

    private static void fillRealDerivativeTable(final double expLambda,
                                                final double[] moments,
                                                final double[] derivatives) {
        for (int leftOrder = 0; leftOrder <= MAX_SMALL_ROOT_DERIVATIVE_ORDER; ++leftOrder) {
            final int derivativeOffset = leftOrder * (MAX_SMALL_ROOT_DERIVATIVE_ORDER + 1);
            for (int rightOrder = 0; rightOrder <= MAX_SMALL_ROOT_DERIVATIVE_ORDER; ++rightOrder) {
                double total = 0.0;
                for (int i = 0; i <= leftOrder; ++i) {
                    final double coefficient = ((i & 1) == 0 ? 1.0 : -1.0) * BINOMIAL[leftOrder][i];
                    total += coefficient * moments[rightOrder + i];
                }
                derivatives[derivativeOffset + rightOrder] = expLambda * total;
            }
        }
    }

    private static void fillRealMomentTable(final double c,
                                            final double[] moments) {
        if (Math.abs(c) < MOMENT_SERIES_CUTOFF) {
            fillRealMomentTableSeries(c, moments);
            return;
        }

        final double expC = Math.exp(c);
        final double expOverC = expC / c;
        moments[0] = (expC - 1.0) / c;
        for (int power = 1; power <= MAX_SMALL_ROOT_MOMENT_POWER; ++power) {
            moments[power] = expOverC - (power * moments[power - 1]) / c;
        }
    }

    private static void fillRealMomentTableSeries(final double c,
                                                  final double[] moments) {
        for (int power = 0; power <= MAX_SMALL_ROOT_MOMENT_POWER; ++power) {
            moments[power] = 0.0;
        }

        double cPower = 1.0;
        double factorial = 1.0;
        for (int n = 0; n < MOMENT_SERIES_TERMS; ++n) {
            if (n > 0) {
                factorial *= n;
                cPower *= c;
            }
            final double base = cPower / factorial;
            for (int power = 0; power <= MAX_SMALL_ROOT_MOMENT_POWER; ++power) {
                moments[power] += base / (n + power + 1.0);
            }
        }
    }

    private static void fillMixedDerivative(final double lambdaRe,
                                            final double lambdaIm,
                                            final double muRe,
                                            final double muIm,
                                            final int leftOrder,
                                            final int rightOrder,
                                            final MutableComplex out,
                                            final ComplexWorkspace workspace) {
        final double cRe = muRe - lambdaRe;
        final double cIm = muIm - lambdaIm;
        double totalRe = 0.0;
        double totalIm = 0.0;
        final MutableComplex moment = workspace.c5;
        for (int i = 0; i <= leftOrder; ++i) {
            final double coefficient = ((i & 1) == 0 ? 1.0 : -1.0) * choose(leftOrder, i);
            fillIntegralMoment(cRe, cIm, rightOrder + i, moment, workspace);
            totalRe += coefficient * moment.re;
            totalIm += coefficient * moment.im;
        }

        final double scale = Math.exp(lambdaRe);
        final double cos = Math.cos(lambdaIm);
        final double sin = Math.sin(lambdaIm);
        out.re = scale * (cos * totalRe - sin * totalIm);
        out.im = scale * (sin * totalRe + cos * totalIm);
    }

    private static void fillPhiIntegral(final double lambdaRe,
                                        final double lambdaIm,
                                        final double muRe,
                                        final double muIm,
                                        final MutableComplex out,
                                        final ComplexWorkspace workspace) {
        fillIntegralMoment(muRe - lambdaRe, muIm - lambdaIm, 0, workspace.c5, workspace);
        final double scale = Math.exp(lambdaRe);
        final double cos = Math.cos(lambdaIm);
        final double sin = Math.sin(lambdaIm);
        out.re = scale * (cos * workspace.c5.re - sin * workspace.c5.im);
        out.im = scale * (sin * workspace.c5.re + cos * workspace.c5.im);
    }

    private static void fillIntegralMoment(final double cRe,
                                           final double cIm,
                                           final int power,
                                           final MutableComplex out,
                                           final ComplexWorkspace workspace) {
        if (Math.hypot(cRe, cIm) < MOMENT_SERIES_CUTOFF) {
            fillIntegrateMomentSeries(cRe, cIm, power, out);
            return;
        }

        final double scale = Math.exp(cRe);
        final double cos = Math.cos(cIm);
        final double sin = Math.sin(cIm);
        final double expRe = scale * cos;
        final double expIm = scale * sin;

        divideComplex(expRe - 1.0, expIm, cRe, cIm, out);
        if (power == 0) {
            return;
        }

        final MutableComplex expOverC = workspace.c6;
        final MutableComplex divided = workspace.c7;
        divideComplex(expRe, expIm, cRe, cIm, expOverC);
        for (int n = 1; n <= power; ++n) {
            divideComplex(n * out.re, n * out.im, cRe, cIm, divided);
            out.re = expOverC.re - divided.re;
            out.im = expOverC.im - divided.im;
        }
    }

    private static void fillIntegrateMomentSeries(final double cRe,
                                                  final double cIm,
                                                  final int power,
                                                  final MutableComplex out) {
        double totalRe = 0.0;
        double totalIm = 0.0;
        double cPowerRe = 1.0;
        double cPowerIm = 0.0;
        double factorial = 1.0;
        for (int n = 0; n < MOMENT_SERIES_TERMS; ++n) {
            if (n > 0) {
                factorial *= n;
                final double nextRe = cPowerRe * cRe - cPowerIm * cIm;
                final double nextIm = cPowerRe * cIm + cPowerIm * cRe;
                cPowerRe = nextRe;
                cPowerIm = nextIm;
            }
            final double termScale = 1.0 / (factorial * (n + power + 1.0));
            totalRe += termScale * cPowerRe;
            totalIm += termScale * cPowerIm;
        }
        out.re = totalRe;
        out.im = totalIm;
    }

    private static boolean hasEqualDiagonals(final double[] block) {
        return Math.abs(block[0] - block[3]) < EPS;
    }

    private static double factorial(final int n) {
        double result = 1.0;
        for (int i = 2; i <= n; ++i) {
            result *= i;
        }
        return result;
    }

    private static double[] buildInverseFactorials(final int maxN) {
        final double[] inverse = new double[maxN + 1];
        double factorial = 1.0;
        inverse[0] = 1.0;
        for (int n = 1; n <= maxN; ++n) {
            factorial *= n;
            inverse[n] = 1.0 / factorial;
        }
        return inverse;
    }

    private static double[][] buildBinomialTable(final int maxN) {
        final double[][] binomial = new double[maxN + 1][maxN + 1];
        for (int n = 0; n <= maxN; ++n) {
            binomial[n][0] = 1.0;
            binomial[n][n] = 1.0;
            for (int k = 1; k < n; ++k) {
                binomial[n][k] = binomial[n - 1][k - 1] + binomial[n - 1][k];
            }
        }
        return binomial;
    }

    private static double choose(final int n,
                                 final int k) {
        if (k < 0 || k > n) {
            return 0.0;
        }
        final int effectiveK = Math.min(k, n - k);
        double result = 1.0;
        for (int i = 1; i <= effectiveK; ++i) {
            result *= (n - effectiveK + i);
            result /= i;
        }
        return result;
    }

    private static void divideComplex(final double numeratorRe,
                                      final double numeratorIm,
                                      final double denominatorRe,
                                      final double denominatorIm,
                                      final MutableComplex out) {
        final double denom = denominatorRe * denominatorRe + denominatorIm * denominatorIm;
        out.re = (numeratorRe * denominatorRe + numeratorIm * denominatorIm) / denom;
        out.im = (numeratorIm * denominatorRe - numeratorRe * denominatorIm) / denom;
    }

    private static void multiplyComplex(final double leftRe,
                                        final double leftIm,
                                        final double rightRe,
                                        final double rightIm,
                                        final MutableComplex out) {
        out.re = leftRe * rightRe - leftIm * rightIm;
        out.im = leftRe * rightIm + leftIm * rightRe;
    }

    private static void fillTwoByTwoBasisMatrix(final int basis,
                                                final double[] out) {
        out[0] = 0.0;
        out[1] = 0.0;
        out[2] = 0.0;
        out[3] = 0.0;
        out[basis] = 1.0;
    }

    private static void multiply2x2(final double[] left,
                                    final double[] right,
                                    final double[] out) {
        final double r00 = right[0];
        final double r01 = right[1];
        final double r10 = right[2];
        final double r11 = right[3];
        out[0] = left[0] * r00 + left[1] * r10;
        out[1] = left[0] * r01 + left[1] * r11;
        out[2] = left[2] * r00 + left[3] * r10;
        out[3] = left[2] * r01 + left[3] * r11;
    }

    private static void fillTwoByTwoRhs(final double[] leftExp,
                                        final double[] rightExp,
                                        final double[] basis,
                                        final double[] out) {
        final double[] leftProduct = new double[4];
        final double[] rightProduct = new double[4];
        multiply2x2(leftExp, basis, leftProduct);
        multiply2x2(basis, rightExp, rightProduct);
        out[0] = leftProduct[0] - rightProduct[0];
        out[1] = leftProduct[1] - rightProduct[1];
        out[2] = leftProduct[2] - rightProduct[2];
        out[3] = leftProduct[3] - rightProduct[3];
    }

    private static void fillTwoByTwoSylvesterOperator(final double[] left,
                                                      final double[] right,
                                                      final double[] out) {
        final double a00 = left[0];
        final double a01 = left[1];
        final double a10 = left[2];
        final double a11 = left[3];

        final double b00 = right[0];
        final double b01 = right[1];
        final double b10 = right[2];
        final double b11 = right[3];

        out[0] = a00 - b00;
        out[1] = -b10;
        out[2] = a01;
        out[3] = 0.0;

        out[4] = -b01;
        out[5] = a00 - b11;
        out[6] = 0.0;
        out[7] = a01;

        out[8] = a10;
        out[9] = 0.0;
        out[10] = a11 - b00;
        out[11] = -b10;

        out[12] = 0.0;
        out[13] = a10;
        out[14] = -b01;
        out[15] = a11 - b11;
    }

    private static void solve4x4(final double[] matrix,
                                 final double[] rhs,
                                 final double[] solution) {
        final double[] augmented = new double[20];
        for (int row = 0; row < 4; ++row) {
            final int matrixOffset = 4 * row;
            final int augmentedOffset = 5 * row;
            augmented[augmentedOffset] = matrix[matrixOffset];
            augmented[augmentedOffset + 1] = matrix[matrixOffset + 1];
            augmented[augmentedOffset + 2] = matrix[matrixOffset + 2];
            augmented[augmentedOffset + 3] = matrix[matrixOffset + 3];
            augmented[augmentedOffset + 4] = rhs[row];
        }

        for (int col = 0; col < 4; ++col) {
            int pivot = col;
            double maxAbs = Math.abs(augmented[5 * col + col]);
            for (int row = col + 1; row < 4; ++row) {
                final double value = Math.abs(augmented[5 * row + col]);
                if (value > maxAbs) {
                    maxAbs = value;
                    pivot = row;
                }
            }

            if (maxAbs < EPS) {
                throw new IllegalStateException("Degenerate 2x2 Sylvester block pair");
            }

            if (pivot != col) {
                swapRows(augmented, col, pivot);
            }

            final int pivotOffset = 5 * col;
            final double pivotValue = augmented[pivotOffset + col];
            for (int j = col; j < 5; ++j) {
                augmented[pivotOffset + j] /= pivotValue;
            }

            for (int row = 0; row < 4; ++row) {
                if (row == col) {
                    continue;
                }
                final int rowOffset = 5 * row;
                final double factor = augmented[rowOffset + col];
                if (factor == 0.0) {
                    continue;
                }
                for (int j = col; j < 5; ++j) {
                    augmented[rowOffset + j] -= factor * augmented[pivotOffset + j];
                }
            }
        }

        solution[0] = augmented[4];
        solution[1] = augmented[9];
        solution[2] = augmented[14];
        solution[3] = augmented[19];
    }

    private static void swapRows(final double[] augmented,
                                 final int firstRow,
                                 final int secondRow) {
        final int firstOffset = 5 * firstRow;
        final int secondOffset = 5 * secondRow;
        for (int j = 0; j < 5; ++j) {
            final double tmp = augmented[firstOffset + j];
            augmented[firstOffset + j] = augmented[secondOffset + j];
            augmented[secondOffset + j] = tmp;
        }
    }

    private static double dividedDifferenceExp(final double left,
                                               final double right) {
        if (Math.abs(left - right) < EPS) {
            return Math.exp(left);
        }
        return (Math.exp(left) - Math.exp(right)) / (left - right);
    }

    private static void fillPhi2x2(final double[] matrix,
                                   final double[] out) {
        final double a = matrix[0];
        final double b = matrix[1];
        final double c = matrix[2];
        final double d = matrix[3];

        final double tau = 0.5 * (a + d);
        final double diff = 0.5 * (a - d);
        final double delta2 = diff * diff + b * c;

        final double i0;
        final double i1;

        if (Math.abs(delta2) < EPS) {
            i0 = scalarExpIntegral(tau);
            i1 = weightedExpIntegral(tau);
        } else if (delta2 > 0.0) {
            final double delta = Math.sqrt(delta2);
            final double plus = scalarExpIntegral(tau + delta);
            final double minus = scalarExpIntegral(tau - delta);
            i0 = 0.5 * (plus + minus);
            i1 = 0.5 * (plus - minus) / delta;
        } else {
            final double omega = Math.sqrt(-delta2);
            final double denom = tau * tau + omega * omega;
            final double expTau = Math.exp(tau);
            final double cos = Math.cos(omega);
            final double sin = Math.sin(omega);
            i0 = (expTau * (tau * cos + omega * sin) - tau) / denom;
            i1 = (expTau * (tau * sin - omega * cos) + omega) / (omega * denom);
        }

        final double g1 = i1;
        final double g0 = i0 - tau * g1;

        out[0] = g0 + g1 * a;
        out[1] = g1 * b;
        out[2] = g1 * c;
        out[3] = g0 + g1 * d;
    }

    private static double scalarExpIntegral(final double value) {
        if (Math.abs(value) < EPS) {
            return 1.0;
        }
        return Math.expm1(value) / value;
    }

    private static double weightedExpIntegral(final double value) {
        if (Math.abs(value) < 1.0e-8) {
            return 0.5 + value / 3.0 + value * value / 8.0;
        }
        return (Math.exp(value) * (value - 1.0) + 1.0) / (value * value);
    }

    private static double[] buildTwoByTwoCoefficientSpectral(final double[] leftBlock,
                                                             final double[] rightBlock,
                                                             final Spectrum2x2 leftSpectrum,
                                                             final Spectrum2x2 rightSpectrum) {
        final ComplexNumber[] leftProjector0 =
                buildProjector(leftBlock, leftSpectrum.lambda[0], leftSpectrum.lambda[1]);
        final ComplexNumber[] leftProjector1 =
                buildProjector(leftBlock, leftSpectrum.lambda[1], leftSpectrum.lambda[0]);
        final ComplexNumber[] rightProjector0 =
                buildProjector(rightBlock, rightSpectrum.lambda[0], rightSpectrum.lambda[1]);
        final ComplexNumber[] rightProjector1 =
                buildProjector(rightBlock, rightSpectrum.lambda[1], rightSpectrum.lambda[0]);

        final ComplexNumber[][] leftProjectors = new ComplexNumber[][]{leftProjector0, leftProjector1};
        final ComplexNumber[][] rightProjectors = new ComplexNumber[][]{rightProjector0, rightProjector1};

        final double[] coefficients = new double[16];
        final double[] basis = new double[4];
        for (int basisIndex = 0; basisIndex < 4; ++basisIndex) {
            fillTwoByTwoBasisMatrix(basisIndex, basis);
            ComplexNumber[] blockSum = zeroComplexMatrix();
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    final ComplexNumber weight =
                            dividedDifferenceExp(leftSpectrum.lambda[i], rightSpectrum.lambda[j]);
                    final ComplexNumber[] term =
                            multiplyComplexMatrix(
                                    multiplyComplexReal(leftProjectors[i], basis),
                                    rightProjectors[j]);
                    addScaledComplexMatrix(blockSum, term, weight);
                }
            }
            coefficients[basisIndex] = blockSum[0].re;
            coefficients[4 + basisIndex] = blockSum[1].re;
            coefficients[8 + basisIndex] = blockSum[2].re;
            coefficients[12 + basisIndex] = blockSum[3].re;
        }
        return coefficients;
    }

    private static double[] buildTwoByTwoCoefficientRepeatedSame(final double[] block) {
        final double lambda = 0.5 * (block[0] + block[3]);
        final double expLambda = Math.exp(lambda);
        final double[] nilpotent = new double[]{
                block[0] - lambda,
                block[1],
                block[2],
                block[3] - lambda
        };

        final double[] coefficients = new double[16];
        final double[] basis = new double[4];
        final double[] leftTerm = new double[4];
        final double[] rightTerm = new double[4];
        final double[] mixedTerm = new double[4];

        for (int basisIndex = 0; basisIndex < 4; ++basisIndex) {
            fillTwoByTwoBasisMatrix(basisIndex, basis);
            multiply2x2(nilpotent, basis, leftTerm);
            multiply2x2(basis, nilpotent, rightTerm);
            multiply2x2(leftTerm, nilpotent, mixedTerm);

            coefficients[basisIndex] =
                    expLambda * (basis[0] + 0.5 * leftTerm[0] + 0.5 * rightTerm[0] + (1.0 / 6.0) * mixedTerm[0]);
            coefficients[4 + basisIndex] =
                    expLambda * (basis[1] + 0.5 * leftTerm[1] + 0.5 * rightTerm[1] + (1.0 / 6.0) * mixedTerm[1]);
            coefficients[8 + basisIndex] =
                    expLambda * (basis[2] + 0.5 * leftTerm[2] + 0.5 * rightTerm[2] + (1.0 / 6.0) * mixedTerm[2]);
            coefficients[12 + basisIndex] =
                    expLambda * (basis[3] + 0.5 * leftTerm[3] + 0.5 * rightTerm[3] + (1.0 / 6.0) * mixedTerm[3]);
        }

        return coefficients;
    }

    private static boolean matricesEqual(final double[] left,
                                         final double[] right) {
        for (int i = 0; i < 4; ++i) {
            if (Math.abs(left[i] - right[i]) > EPS) {
                return false;
            }
        }
        return true;
    }

    private static ComplexNumber[] buildProjector(final double[] block,
                                                  final ComplexNumber lambda,
                                                  final ComplexNumber other) {
        final ComplexNumber denominator = lambda.subtract(other);
        return new ComplexNumber[]{
                ComplexNumber.real(block[0]).subtract(other).divide(denominator),
                ComplexNumber.real(block[1]).divide(denominator),
                ComplexNumber.real(block[2]).divide(denominator),
                ComplexNumber.real(block[3]).subtract(other).divide(denominator)
        };
    }

    private static ComplexNumber[] multiplyComplexReal(final ComplexNumber[] left,
                                                       final double[] right) {
        return new ComplexNumber[]{
                left[0].scale(right[0]).add(left[1].scale(right[2])),
                left[0].scale(right[1]).add(left[1].scale(right[3])),
                left[2].scale(right[0]).add(left[3].scale(right[2])),
                left[2].scale(right[1]).add(left[3].scale(right[3]))
        };
    }

    private static ComplexNumber[] multiplyComplexMatrix(final ComplexNumber[] left,
                                                         final ComplexNumber[] right) {
        return new ComplexNumber[]{
                left[0].multiply(right[0]).add(left[1].multiply(right[2])),
                left[0].multiply(right[1]).add(left[1].multiply(right[3])),
                left[2].multiply(right[0]).add(left[3].multiply(right[2])),
                left[2].multiply(right[1]).add(left[3].multiply(right[3]))
        };
    }

    private static ComplexNumber[] zeroComplexMatrix() {
        return new ComplexNumber[]{
                ComplexNumber.ZERO,
                ComplexNumber.ZERO,
                ComplexNumber.ZERO,
                ComplexNumber.ZERO
        };
    }

    private static void addScaledComplexMatrix(final ComplexNumber[] accumulator,
                                               final ComplexNumber[] term,
                                               final ComplexNumber scale) {
        accumulator[0] = accumulator[0].add(scale.multiply(term[0]));
        accumulator[1] = accumulator[1].add(scale.multiply(term[1]));
        accumulator[2] = accumulator[2].add(scale.multiply(term[2]));
        accumulator[3] = accumulator[3].add(scale.multiply(term[3]));
    }

    private static final class Spectrum2x2 {
        private final boolean distinct;
        private final ComplexNumber[] lambda;

        private Spectrum2x2(final boolean distinct,
                            final ComplexNumber[] lambda) {
            this.distinct = distinct;
            this.lambda = lambda;
        }

        private static Spectrum2x2 decompose(final double[] block) {
            final double trace = block[0] + block[3];
            final double determinant = block[0] * block[3] - block[1] * block[2];
            final double discriminant = trace * trace - 4.0 * determinant;
            if (Math.abs(discriminant) < EPS) {
                final ComplexNumber repeated = ComplexNumber.real(0.5 * trace);
                return new Spectrum2x2(false, new ComplexNumber[]{repeated, repeated});
            }
            if (discriminant > 0.0) {
                final double root = Math.sqrt(discriminant);
                return new Spectrum2x2(
                        true,
                        new ComplexNumber[]{
                                ComplexNumber.real(0.5 * (trace + root)),
                                ComplexNumber.real(0.5 * (trace - root))
                        });
            }
            final double imag = 0.5 * Math.sqrt(-discriminant);
            final double real = 0.5 * trace;
            return new Spectrum2x2(
                    true,
                    new ComplexNumber[]{
                            new ComplexNumber(real, imag),
                            new ComplexNumber(real, -imag)
                    });
        }
    }

    private static final class EqualDiagonalCoefficients {
        private final double alpha;
        private final double beta;
        private final double gamma;
        private final double eta;

        private EqualDiagonalCoefficients(final double alpha,
                                          final double beta,
                                          final double gamma,
                                          final double eta) {
            this.alpha = alpha;
            this.beta = beta;
            this.gamma = gamma;
            this.eta = eta;
        }
    }

    private static final class MutableComplex {
        private double re;
        private double im;
    }

    private static final class ComplexWorkspace {
        private final MutableComplex c0 = new MutableComplex();
        private final MutableComplex c1 = new MutableComplex();
        private final MutableComplex c2 = new MutableComplex();
        private final MutableComplex c3 = new MutableComplex();
        private final MutableComplex c4 = new MutableComplex();
        private final MutableComplex c5 = new MutableComplex();
        private final MutableComplex c6 = new MutableComplex();
        private final MutableComplex c7 = new MutableComplex();
    }

    private static final class RealSeriesWorkspace {
        private final double[] moments = new double[MAX_SMALL_ROOT_MOMENT_POWER + 1];
        private final double[] derivatives =
                new double[(MAX_SMALL_ROOT_DERIVATIVE_ORDER + 1) * (MAX_SMALL_ROOT_DERIVATIVE_ORDER + 1)];
    }

    private static final class ComplexNumber {
        private static final ComplexNumber ZERO = new ComplexNumber(0.0, 0.0);
        private static final ComplexNumber ONE = new ComplexNumber(1.0, 0.0);

        private final double re;
        private final double im;

        private ComplexNumber(final double re,
                              final double im) {
            this.re = re;
            this.im = im;
        }

        private static ComplexNumber real(final double value) {
            return new ComplexNumber(value, 0.0);
        }

        private static ComplexNumber sqrtOfReal(final double value) {
            if (value >= 0.0) {
                return new ComplexNumber(Math.sqrt(value), 0.0);
            }
            return new ComplexNumber(0.0, Math.sqrt(-value));
        }

        private ComplexNumber add(final ComplexNumber other) {
            return new ComplexNumber(re + other.re, im + other.im);
        }

        private ComplexNumber subtract(final ComplexNumber other) {
            return new ComplexNumber(re - other.re, im - other.im);
        }

        private ComplexNumber scale(final double factor) {
            return new ComplexNumber(re * factor, im * factor);
        }

        private ComplexNumber pow(final int exponent) {
            if (exponent < 0) {
                throw new IllegalArgumentException("Negative exponent: " + exponent);
            }
            ComplexNumber result = ONE;
            for (int i = 0; i < exponent; ++i) {
                result = result.multiply(this);
            }
            return result;
        }

        private ComplexNumber multiply(final ComplexNumber other) {
            return new ComplexNumber(
                    re * other.re - im * other.im,
                    re * other.im + im * other.re);
        }

        private ComplexNumber divide(final ComplexNumber other) {
            final double denom = other.re * other.re + other.im * other.im;
            return new ComplexNumber(
                    (re * other.re + im * other.im) / denom,
                    (im * other.re - re * other.im) / denom);
        }

        private double absSquared() {
            return re * re + im * im;
        }

        private double abs() {
            return Math.sqrt(absSquared());
        }

        private ComplexNumber exp() {
            final double scale = Math.exp(re);
            return new ComplexNumber(scale * Math.cos(im), scale * Math.sin(im));
        }
    }

    private static ComplexNumber dividedDifferenceExp(final ComplexNumber left,
                                                      final ComplexNumber right) {
        final ComplexNumber delta = left.subtract(right);
        if (delta.absSquared() < EPS * EPS) {
            return left.exp();
        }
        return left.exp().subtract(right.exp()).divide(delta);
    }
}
