package dr.evomodel.treedatalikelihood.continuous.backprop;

import org.ejml.data.DenseMatrix64F;

import java.util.concurrent.atomic.AtomicLong;

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
    private static final ThreadLocal<BuildWorkspace> BUILD_WORKSPACE =
            ThreadLocal.withInitial(BuildWorkspace::new);
    private static final ThreadLocal<ComplexWorkspace> COMPLEX_WORKSPACE =
            ThreadLocal.withInitial(ComplexWorkspace::new);
    private static final ThreadLocal<RealSeriesWorkspace> REAL_SERIES_WORKSPACE =
            ThreadLocal.withInitial(RealSeriesWorkspace::new);
    private static final AtomicLong PARAMETER_UPDATES = new AtomicLong();
    private static final AtomicLong TIME_EVALUATIONS = new AtomicLong();
    private static final AtomicLong EQUAL_DIAGONAL_EVALUATIONS = new AtomicLong();
    private static final AtomicLong GENERIC_SOLVE4X4_CALLS = new AtomicLong();

    private BlockDiagonalFrechetExactPlan() {
    }

    static Plan createPlan(final BlockDiagonalExpSolver.BlockStructure structure) {
        return new Plan(structure);
    }

    static void resetInstrumentation() {
        PARAMETER_UPDATES.set(0L);
        TIME_EVALUATIONS.set(0L);
        EQUAL_DIAGONAL_EVALUATIONS.set(0L);
        GENERIC_SOLVE4X4_CALLS.set(0L);
    }

    static long getParameterUpdateCount() {
        return PARAMETER_UPDATES.get();
    }

    static long getTimeEvaluationCount() {
        return TIME_EVALUATIONS.get();
    }

    static long getEqualDiagonalEvaluationCount() {
        return EQUAL_DIAGONAL_EVALUATIONS.get();
    }

    static long getGenericSolve4x4CallCount() {
        return GENERIC_SOLVE4X4_CALLS.get();
    }

    private static Entry[] buildEntries(final BlockDiagonalExpSolver.BlockStructure structure) {
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
                        rightDim);
            }
        }
        return entries;
    }

    static final class Plan {
        private final Entry[] entries;
        private final int dimension;

        private Plan(final BlockDiagonalExpSolver.BlockStructure structure) {
            this.dimension = structure.getDim();
            this.entries = buildEntries(structure);
        }

        void updateParameters(final double[] blockDParams) {
            PARAMETER_UPDATES.incrementAndGet();
            for (final Entry entry : entries) {
                entry.kernel.updateParameters(blockDParams, dimension);
            }
        }

        void evaluate(final double t) {
            TIME_EVALUATIONS.incrementAndGet();
            final BuildWorkspace buildWorkspace = BUILD_WORKSPACE.get();
            for (final Entry entry : entries) {
                entry.kernel.fill(t, entry.coefficients, buildWorkspace);
            }
        }

        void rebuild(final double[] blockDParams,
                     final double t) {
            updateParameters(blockDParams);
            evaluate(t);
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
        private final EntryKernel kernel;

        private Entry(final int leftStart,
                      final int leftDim,
                      final int rightStart,
                      final int rightDim) {
            this.leftStart = leftStart;
            this.leftDim = leftDim;
            this.rightStart = rightStart;
            this.rightDim = rightDim;
            this.size = leftDim * rightDim;
            this.coefficients = new double[size * size];
            this.kernel = createKernel(leftStart, leftDim, rightStart, rightDim);
        }
    }

    private interface EntryKernel {
        void updateParameters(double[] blockDParams, int dimension);

        void fill(double t, double[] out, BuildWorkspace workspace);
    }

    private static EntryKernel createKernel(final int leftStart,
                                            final int leftDim,
                                            final int rightStart,
                                            final int rightDim) {
        if (leftDim == 2 && rightDim == 2) {
            return new TwoByTwoBlockPairKernel(leftStart, rightStart);
        }
        return new ExistingCoefficientKernel(leftStart, leftDim, rightStart, rightDim);
    }

    private static final class ExistingCoefficientKernel implements EntryKernel {
        private final int leftStart;
        private final int leftDim;
        private final int rightStart;
        private final int rightDim;
        private double[] blockDParams;
        private int dimension;

        private ExistingCoefficientKernel(final int leftStart,
                                          final int leftDim,
                                          final int rightStart,
                                          final int rightDim) {
            this.leftStart = leftStart;
            this.leftDim = leftDim;
            this.rightStart = rightStart;
            this.rightDim = rightDim;
        }

        @Override
        public void updateParameters(final double[] blockDParams,
                                     final int dimension) {
            this.blockDParams = blockDParams;
            this.dimension = dimension;
        }

        @Override
        public void fill(final double t,
                         final double[] out,
                         final BuildWorkspace workspace) {
            fillCoefficients(
                    blockDParams,
                    dimension,
                    leftStart,
                    leftDim,
                    rightStart,
                    rightDim,
                    t,
                    out,
                    workspace);
        }
    }

    private static final class TwoByTwoBlockPairKernel implements EntryKernel {
        private final int leftStart;
        private final int rightStart;
        private final EqualDiagonalBlockBlockKernel equalDiagonalKernel =
                new EqualDiagonalBlockBlockKernel();
        private final double[] leftBlock = new double[4];
        private final double[] rightBlock = new double[4];
        private final double[] leftExp = new double[4];
        private final double[] rightExp = new double[4];
        private double leftDiagonal;
        private double leftUpper;
        private double leftLower;
        private double leftOtherDiagonal;
        private double rightDiagonal;
        private double rightUpper;
        private double rightLower;
        private double rightOtherDiagonal;
        private double leftProduct;
        private double rightProduct;
        private boolean equalDiagonalPair;

        private TwoByTwoBlockPairKernel(final int leftStart,
                                        final int rightStart) {
            this.leftStart = leftStart;
            this.rightStart = rightStart;
        }

        @Override
        public void updateParameters(final double[] blockDParams,
                                     final int dimension) {
            final int upperOffset = dimension;
            final int lowerOffset = dimension + (dimension - 1);
            leftDiagonal = blockDParams[leftStart];
            leftUpper = blockDParams[upperOffset + leftStart];
            leftLower = blockDParams[lowerOffset + leftStart];
            leftOtherDiagonal = blockDParams[leftStart + 1];
            rightDiagonal = blockDParams[rightStart];
            rightUpper = blockDParams[upperOffset + rightStart];
            rightLower = blockDParams[lowerOffset + rightStart];
            rightOtherDiagonal = blockDParams[rightStart + 1];
            leftProduct = leftUpper * leftLower;
            rightProduct = rightUpper * rightLower;
            equalDiagonalPair =
                    Math.abs(leftDiagonal - leftOtherDiagonal) < EPS
                            && Math.abs(rightDiagonal - rightOtherDiagonal) < EPS;
            if (equalDiagonalPair) {
                equalDiagonalKernel.update(
                        leftDiagonal,
                        leftUpper,
                        leftLower,
                        rightDiagonal,
                        rightUpper,
                        rightLower,
                        leftProduct,
                        rightProduct);
            }
        }

        @Override
        public void fill(final double t,
                         final double[] out,
                         final BuildWorkspace workspace) {
            if (equalDiagonalPair) {
                fillEqualDiagonalCoefficient(t, out, workspace);
                return;
            }

            fillScaledBlocks(t);
            fillSmallExponential(leftBlock, 2, leftExp);
            fillSmallExponential(rightBlock, 2, rightExp);
            fillTwoByTwoCoefficientGeneric(leftBlock, rightBlock, leftExp, rightExp, out, workspace);
        }

        private void fillEqualDiagonalCoefficient(final double t,
                                                  final double[] out,
                                                  final BuildWorkspace workspace) {
            EQUAL_DIAGONAL_EVALUATIONS.incrementAndGet();
            equalDiagonalKernel.fill(t, out, workspace);
        }

        private void fillScaledBlocks(final double t) {
            final double scale = -t;
            leftBlock[0] = scale * leftDiagonal;
            leftBlock[1] = scale * leftUpper;
            leftBlock[2] = scale * leftLower;
            leftBlock[3] = scale * leftOtherDiagonal;
            rightBlock[0] = scale * rightDiagonal;
            rightBlock[1] = scale * rightUpper;
            rightBlock[2] = scale * rightLower;
            rightBlock[3] = scale * rightOtherDiagonal;
        }
    }

    private static final class EqualDiagonalBlockBlockKernel {
        private double leftDiagonal;
        private double leftUpper;
        private double leftLower;
        private double rightDiagonal;
        private double rightUpper;
        private double rightLower;
        private double leftProduct;
        private double rightProduct;

        void update(final double leftDiagonal,
                    final double leftUpper,
                    final double leftLower,
                    final double rightDiagonal,
                    final double rightUpper,
                    final double rightLower,
                    final double leftProduct,
                    final double rightProduct) {
            this.leftDiagonal = leftDiagonal;
            this.leftUpper = leftUpper;
            this.leftLower = leftLower;
            this.rightDiagonal = rightDiagonal;
            this.rightUpper = rightUpper;
            this.rightLower = rightLower;
            this.leftProduct = leftProduct;
            this.rightProduct = rightProduct;
        }

        void fill(final double t,
                  final double[] out,
                  final BuildWorkspace workspace) {
            final double scale = -t;
            final double a = scale * leftDiagonal;
            final double b = scale * rightDiagonal;
            final double scaledLeftUpper = scale * leftUpper;
            final double scaledLeftLower = scale * leftLower;
            final double scaledRightUpper = scale * rightUpper;
            final double scaledRightLower = scale * rightLower;
            final double t2 = t * t;
            final double scaledLeftProduct = t2 * leftProduct;
            final double scaledRightProduct = t2 * rightProduct;
            final double rRe = scaledLeftProduct >= 0.0 ? Math.sqrt(scaledLeftProduct) : 0.0;
            final double rIm = scaledLeftProduct >= 0.0 ? 0.0 : Math.sqrt(-scaledLeftProduct);
            final double qRe = scaledRightProduct >= 0.0 ? Math.sqrt(scaledRightProduct) : 0.0;
            final double qIm = scaledRightProduct >= 0.0 ? 0.0 : Math.sqrt(-scaledRightProduct);
            final double rAbs = Math.sqrt(Math.abs(scaledLeftProduct));
            final double qAbs = Math.sqrt(Math.abs(scaledRightProduct));
            final ComplexWorkspace complexWorkspace = COMPLEX_WORKSPACE.get();
            final PreparedCoefficients prepared = workspace.prepared;

            if (rAbs < SMALL_ROOT_SERIES_CUTOFF && qAbs < SMALL_ROOT_SERIES_CUTOFF) {
                fillCoefficientsBivariateSmallRootSeries(a, b, scaledLeftProduct, scaledRightProduct, prepared);
            } else if (rAbs < SMALL_ROOT_SERIES_CUTOFF) {
                fillCoefficientsSmallLeftRootSeries(
                        a, b, qRe, qIm, scaledLeftProduct, prepared, complexWorkspace);
            } else if (qAbs < SMALL_ROOT_SERIES_CUTOFF) {
                fillCoefficientsSmallRightRootSeries(
                        a, b, rRe, rIm, scaledRightProduct, prepared, complexWorkspace);
            } else {
                fillCoefficientsDistinct(a, b, rRe, rIm, qRe, qIm, prepared, complexWorkspace);
            }

            fillTwoByTwoCoefficientFromPrepared(
                    prepared,
                    scaledLeftUpper,
                    scaledLeftLower,
                    scaledRightUpper,
                    scaledRightLower,
                    out,
                    workspace);
        }
    }

    private static void fillCoefficients(final double[] blockDParams,
                                         final int dimension,
                                         final int leftStart,
                                         final int leftDim,
                                         final int rightStart,
                                         final int rightDim,
                                         final double t,
                                         final double[] out,
                                         final BuildWorkspace workspace) {
        final double[] left = workspace.left;
        final double[] right = workspace.right;
        final double[] leftExp = workspace.leftExp;
        final double[] rightExp = workspace.rightExp;

        fillScaledBlock(blockDParams, dimension, leftStart, leftDim, -t, left);
        fillScaledBlock(blockDParams, dimension, rightStart, rightDim, -t, right);
        fillSmallExponential(left, leftDim, leftExp);
        fillSmallExponential(right, rightDim, rightExp);

        if (leftDim == 1 && rightDim == 1) {
            out[0] = dividedDifferenceExp(left[0], right[0]);
            return;
        }
        if (leftDim == 1 && rightDim == 2) {
            fillOneByTwoCoefficient(left[0], right, out, workspace);
            return;
        }
        if (leftDim == 2 && rightDim == 1) {
            fillTwoByOneCoefficient(left, right[0], out, workspace);
            return;
        }
        if (leftDim == 2 && rightDim == 2) {
            fillTwoByTwoCoefficient(left, right, leftExp, rightExp, out, workspace);
            return;
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

    private static void fillOneByTwoCoefficient(final double leftScalar,
                                                final double[] rightBlock,
                                                final double[] out,
                                                final BuildWorkspace workspace) {
        final double[] shifted = workspace.shifted;
        final double[] phi = workspace.phi;
        shifted[0] = rightBlock[0] - leftScalar;
        shifted[1] = rightBlock[1];
        shifted[2] = rightBlock[2];
        shifted[3] = rightBlock[3] - leftScalar;
        fillPhi2x2(shifted, phi);

        final double scale = Math.exp(leftScalar);
        phi[0] *= scale;
        phi[1] *= scale;
        phi[2] *= scale;
        phi[3] *= scale;

        out[0] = phi[0];
        out[1] = phi[2];
        out[2] = phi[1];
        out[3] = phi[3];
    }

    private static void fillTwoByOneCoefficient(final double[] leftBlock,
                                                final double rightScalar,
                                                final double[] out,
                                                final BuildWorkspace workspace) {
        final double[] shifted = workspace.shifted;
        final double[] phi = workspace.phi;
        shifted[0] = leftBlock[0] - rightScalar;
        shifted[1] = leftBlock[1];
        shifted[2] = leftBlock[2];
        shifted[3] = leftBlock[3] - rightScalar;
        fillPhi2x2(shifted, phi);

        final double scale = Math.exp(rightScalar);
        for (int i = 0; i < 4; ++i) {
            out[i] = phi[i] * scale;
        }
    }

    private static void fillTwoByTwoCoefficient(final double[] leftBlock,
                                                final double[] rightBlock,
                                                final double[] leftExp,
                                                final double[] rightExp,
                                                final double[] out,
                                                final BuildWorkspace workspace) {
        if (hasEqualDiagonals(leftBlock) && hasEqualDiagonals(rightBlock)) {
            fillTwoByTwoCoefficientEqualDiagonal(leftBlock, rightBlock, out, workspace);
            return;
        }
        fillTwoByTwoCoefficientGeneric(leftBlock, rightBlock, leftExp, rightExp, out, workspace);
    }

    private static void fillTwoByTwoCoefficientEqualDiagonal(final double[] leftBlock,
                                                             final double[] rightBlock,
                                                             final double[] out,
                                                             final BuildWorkspace workspace) {
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
        final ComplexWorkspace complexWorkspace = COMPLEX_WORKSPACE.get();
        final PreparedCoefficients prepared = workspace.prepared;

        if (rAbs < SMALL_ROOT_SERIES_CUTOFF && qAbs < SMALL_ROOT_SERIES_CUTOFF) {
            fillCoefficientsBivariateSmallRootSeries(a, b, leftProduct, rightProduct, prepared);
        } else if (rAbs < SMALL_ROOT_SERIES_CUTOFF) {
            fillCoefficientsSmallLeftRootSeries(a, b, qRe, qIm, leftProduct, prepared, complexWorkspace);
        } else if (qAbs < SMALL_ROOT_SERIES_CUTOFF) {
            fillCoefficientsSmallRightRootSeries(a, b, rRe, rIm, rightProduct, prepared, complexWorkspace);
        } else {
            fillCoefficientsDistinct(a, b, rRe, rIm, qRe, qIm, prepared, complexWorkspace);
        }

        fillTwoByTwoCoefficientFromPrepared(prepared, leftUpper, leftLower, rightUpper, rightLower, out, workspace);
    }

    private static void fillTwoByTwoCoefficientGeneric(final double[] leftBlock,
                                                       final double[] rightBlock,
                                                       final double[] leftExp,
                                                       final double[] rightExp,
                                                       final double[] out,
                                                       final BuildWorkspace workspace) {
        try {
            final double[] sylvester = workspace.sylvester;
            fillTwoByTwoSylvesterOperator(leftBlock, rightBlock, sylvester);

            final double[] rhs = workspace.rhs;
            final double[] solution = workspace.solution;
            final double[] basisMatrix = workspace.basis;
            final double[] augmented = workspace.augmented;

            for (int basis = 0; basis < 4; ++basis) {
                fillTwoByTwoBasisMatrix(basis, basisMatrix);
                fillTwoByTwoRhs(
                        leftExp,
                        rightExp,
                        basisMatrix,
                        rhs,
                        workspace.leftProduct,
                        workspace.rightProduct);
                solve4x4(sylvester, rhs, solution, augmented);
                out[basis] = solution[0];
                out[4 + basis] = solution[1];
                out[8 + basis] = solution[2];
                out[12 + basis] = solution[3];
            }
        } catch (final IllegalStateException degenerate) {
            fillSpectrum2x2(leftBlock, workspace.leftSpectrum);
            fillSpectrum2x2(rightBlock, workspace.rightSpectrum);
            if (workspace.leftSpectrum.distinct && workspace.rightSpectrum.distinct) {
                fillTwoByTwoCoefficientSpectral(leftBlock, rightBlock, out, workspace);
                return;
            }
            if (matricesEqual(leftBlock, rightBlock)) {
                fillTwoByTwoCoefficientRepeatedSame(leftBlock, out, workspace);
                return;
            }
            throw new IllegalStateException(
                    "Degenerate 2x2 Sylvester block pair left=[" + leftBlock[0] + ", " + leftBlock[1] + ", "
                            + leftBlock[2] + ", " + leftBlock[3] + "] right=["
                            + rightBlock[0] + ", " + rightBlock[1] + ", "
                            + rightBlock[2] + ", " + rightBlock[3] + "] leftDistinct="
                            + workspace.leftSpectrum.distinct + " rightDistinct=" + workspace.rightSpectrum.distinct
                            + " equal=" + matricesEqual(leftBlock, rightBlock),
                    degenerate);
        }
    }

    private static void fillTwoByTwoCoefficientFromPrepared(final PreparedCoefficients prepared,
                                                            final double leftUpper,
                                                            final double leftLower,
                                                            final double rightUpper,
                                                            final double rightLower,
                                                            final double[] out,
                                                            final BuildWorkspace workspace) {
        final double[] basis = workspace.basis;
        for (int basisIndex = 0; basisIndex < 4; ++basisIndex) {
            fillTwoByTwoBasisMatrix(basisIndex, basis);
            out[basisIndex] =
                    prepared.alpha * basis[0]
                            + prepared.beta * (leftUpper * basis[2])
                            + prepared.gamma * (basis[1] * rightLower)
                            + prepared.eta * (leftUpper * basis[3] * rightLower);
            out[4 + basisIndex] =
                    prepared.alpha * basis[1]
                            + prepared.beta * (leftUpper * basis[3])
                            + prepared.gamma * (basis[0] * rightUpper)
                            + prepared.eta * (leftUpper * basis[2] * rightUpper);
            out[8 + basisIndex] =
                    prepared.alpha * basis[2]
                            + prepared.beta * (leftLower * basis[0])
                            + prepared.gamma * (basis[3] * rightLower)
                            + prepared.eta * (leftLower * basis[1] * rightLower);
            out[12 + basisIndex] =
                    prepared.alpha * basis[3]
                            + prepared.beta * (leftLower * basis[1])
                            + prepared.gamma * (basis[2] * rightUpper)
                            + prepared.eta * (leftLower * basis[0] * rightUpper);
        }
    }

    private static void fillCoefficientsDistinct(final double a,
                                                 final double b,
                                                 final double rRe,
                                                 final double rIm,
                                                 final double qRe,
                                                 final double qIm,
                                                 final PreparedCoefficients out,
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
        out.alpha = alpha;
        out.beta = beta;
        out.gamma = gamma;
        out.eta = 0.25 * scratch.re;
    }

    private static void fillCoefficientsBivariateSmallRootSeries(final double a,
                                                                 final double b,
                                                                 final double rSquared,
                                                                 final double qSquared,
                                                                 final PreparedCoefficients out) {
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

        out.alpha = alpha;
        out.beta = beta;
        out.gamma = gamma;
        out.eta = eta;
    }

    private static void fillCoefficientsSmallLeftRootSeries(final double a,
                                                            final double b,
                                                            final double qRe,
                                                            final double qIm,
                                                            final double rSquared,
                                                            final PreparedCoefficients out,
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

        out.alpha = alpha;
        out.beta = beta;
        out.gamma = gamma;
        out.eta = eta;
    }

    private static void fillCoefficientsSmallRightRootSeries(final double a,
                                                             final double b,
                                                             final double rRe,
                                                             final double rIm,
                                                             final double qSquared,
                                                             final PreparedCoefficients out,
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

        out.alpha = alpha;
        out.beta = beta;
        out.gamma = gamma;
        out.eta = eta;
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
                                        final double[] out,
                                        final double[] leftProduct,
                                        final double[] rightProduct) {
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
                                 final double[] solution,
                                 final double[] augmented) {
        GENERIC_SOLVE4X4_CALLS.incrementAndGet();
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

    private static void fillTwoByTwoCoefficientSpectral(final double[] leftBlock,
                                                        final double[] rightBlock,
                                                        final double[] out,
                                                        final BuildWorkspace workspace) {
        final MutableSpectrum leftSpectrum = workspace.leftSpectrum;
        final MutableSpectrum rightSpectrum = workspace.rightSpectrum;
        fillProjector(leftBlock, leftSpectrum.lambda0, leftSpectrum.lambda1, workspace.leftProjector0);
        fillProjector(leftBlock, leftSpectrum.lambda1, leftSpectrum.lambda0, workspace.leftProjector1);
        fillProjector(rightBlock, rightSpectrum.lambda0, rightSpectrum.lambda1, workspace.rightProjector0);
        fillProjector(rightBlock, rightSpectrum.lambda1, rightSpectrum.lambda0, workspace.rightProjector1);

        for (int basisIndex = 0; basisIndex < 4; ++basisIndex) {
            fillTwoByTwoBasisMatrix(basisIndex, workspace.basis);
            clearComplexMatrix(workspace.complexAccumulator);
            for (int i = 0; i < 2; ++i) {
                final ComplexMatrix leftProjector = i == 0 ? workspace.leftProjector0 : workspace.leftProjector1;
                final MutableComplex leftLambda = i == 0 ? leftSpectrum.lambda0 : leftSpectrum.lambda1;
                multiplyComplexReal(leftProjector, workspace.basis, workspace.leftTimesBasis);
                for (int j = 0; j < 2; ++j) {
                    final ComplexMatrix rightProjector = j == 0 ? workspace.rightProjector0 : workspace.rightProjector1;
                    final MutableComplex rightLambda = j == 0 ? rightSpectrum.lambda0 : rightSpectrum.lambda1;
                    fillDividedDifferenceExpComplex(
                            leftLambda.re, leftLambda.im,
                            rightLambda.re, rightLambda.im,
                            workspace.complexScalar);
                    multiplyComplexMatrix(workspace.leftTimesBasis, rightProjector, workspace.complexTerm);
                    addScaledComplexMatrix(
                            workspace.complexAccumulator,
                            workspace.complexTerm,
                            workspace.complexScalar.re,
                            workspace.complexScalar.im);
                }
            }
            out[basisIndex] = workspace.complexAccumulator.data[0].re;
            out[4 + basisIndex] = workspace.complexAccumulator.data[1].re;
            out[8 + basisIndex] = workspace.complexAccumulator.data[2].re;
            out[12 + basisIndex] = workspace.complexAccumulator.data[3].re;
        }
    }

    private static void fillTwoByTwoCoefficientRepeatedSame(final double[] block,
                                                            final double[] out,
                                                            final BuildWorkspace workspace) {
        final double lambda = 0.5 * (block[0] + block[3]);
        final double expLambda = Math.exp(lambda);
        final double[] nilpotent = workspace.nilpotent;
        nilpotent[0] = block[0] - lambda;
        nilpotent[1] = block[1];
        nilpotent[2] = block[2];
        nilpotent[3] = block[3] - lambda;

        for (int basisIndex = 0; basisIndex < 4; ++basisIndex) {
            fillTwoByTwoBasisMatrix(basisIndex, workspace.basis);
            multiply2x2(nilpotent, workspace.basis, workspace.leftTerm);
            multiply2x2(workspace.basis, nilpotent, workspace.rightTerm);
            multiply2x2(workspace.leftTerm, nilpotent, workspace.mixedTerm);

            out[basisIndex] =
                    expLambda * (workspace.basis[0]
                    + 0.5 * workspace.leftTerm[0]
                    + 0.5 * workspace.rightTerm[0]
                    + (1.0 / 6.0) * workspace.mixedTerm[0]);
            out[4 + basisIndex] =
                    expLambda * (workspace.basis[1]
                    + 0.5 * workspace.leftTerm[1]
                    + 0.5 * workspace.rightTerm[1]
                    + (1.0 / 6.0) * workspace.mixedTerm[1]);
            out[8 + basisIndex] =
                    expLambda * (workspace.basis[2]
                    + 0.5 * workspace.leftTerm[2]
                    + 0.5 * workspace.rightTerm[2]
                    + (1.0 / 6.0) * workspace.mixedTerm[2]);
            out[12 + basisIndex] =
                    expLambda * (workspace.basis[3]
                    + 0.5 * workspace.leftTerm[3]
                    + 0.5 * workspace.rightTerm[3]
                    + (1.0 / 6.0) * workspace.mixedTerm[3]);
        }
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

    private static void fillSpectrum2x2(final double[] block,
                                        final MutableSpectrum out) {
        final double trace = block[0] + block[3];
        final double determinant = block[0] * block[3] - block[1] * block[2];
        final double discriminant = trace * trace - 4.0 * determinant;
        if (Math.abs(discriminant) < EPS) {
            out.distinct = false;
            setComplex(out.lambda0, 0.5 * trace, 0.0);
            setComplex(out.lambda1, 0.5 * trace, 0.0);
            return;
        }
        out.distinct = true;
        if (discriminant > 0.0) {
            final double root = Math.sqrt(discriminant);
            setComplex(out.lambda0, 0.5 * (trace + root), 0.0);
            setComplex(out.lambda1, 0.5 * (trace - root), 0.0);
            return;
        }
        final double imag = 0.5 * Math.sqrt(-discriminant);
        final double real = 0.5 * trace;
        setComplex(out.lambda0, real, imag);
        setComplex(out.lambda1, real, -imag);
    }

    private static void fillProjector(final double[] block,
                                      final MutableComplex lambda,
                                      final MutableComplex other,
                                      final ComplexMatrix out) {
        final double denominatorRe = lambda.re - other.re;
        final double denominatorIm = lambda.im - other.im;
        divideComplex(block[0] - other.re, -other.im, denominatorRe, denominatorIm, out.data[0]);
        divideComplex(block[1], 0.0, denominatorRe, denominatorIm, out.data[1]);
        divideComplex(block[2], 0.0, denominatorRe, denominatorIm, out.data[2]);
        divideComplex(block[3] - other.re, -other.im, denominatorRe, denominatorIm, out.data[3]);
    }

    private static void multiplyComplexReal(final ComplexMatrix left,
                                            final double[] right,
                                            final ComplexMatrix out) {
        final double r00 = right[0];
        final double r01 = right[1];
        final double r10 = right[2];
        final double r11 = right[3];
        setComplex(
                out.data[0],
                left.data[0].re * r00 + left.data[1].re * r10,
                left.data[0].im * r00 + left.data[1].im * r10);
        setComplex(
                out.data[1],
                left.data[0].re * r01 + left.data[1].re * r11,
                left.data[0].im * r01 + left.data[1].im * r11);
        setComplex(
                out.data[2],
                left.data[2].re * r00 + left.data[3].re * r10,
                left.data[2].im * r00 + left.data[3].im * r10);
        setComplex(
                out.data[3],
                left.data[2].re * r01 + left.data[3].re * r11,
                left.data[2].im * r01 + left.data[3].im * r11);
    }

    private static void multiplyComplexMatrix(final ComplexMatrix left,
                                              final ComplexMatrix right,
                                              final ComplexMatrix out) {
        setComplex(
                out.data[0],
                left.data[0].re * right.data[0].re - left.data[0].im * right.data[0].im
                        + left.data[1].re * right.data[2].re - left.data[1].im * right.data[2].im,
                left.data[0].re * right.data[0].im + left.data[0].im * right.data[0].re
                        + left.data[1].re * right.data[2].im + left.data[1].im * right.data[2].re);
        setComplex(
                out.data[1],
                left.data[0].re * right.data[1].re - left.data[0].im * right.data[1].im
                        + left.data[1].re * right.data[3].re - left.data[1].im * right.data[3].im,
                left.data[0].re * right.data[1].im + left.data[0].im * right.data[1].re
                        + left.data[1].re * right.data[3].im + left.data[1].im * right.data[3].re);
        setComplex(
                out.data[2],
                left.data[2].re * right.data[0].re - left.data[2].im * right.data[0].im
                        + left.data[3].re * right.data[2].re - left.data[3].im * right.data[2].im,
                left.data[2].re * right.data[0].im + left.data[2].im * right.data[0].re
                        + left.data[3].re * right.data[2].im + left.data[3].im * right.data[2].re);
        setComplex(
                out.data[3],
                left.data[2].re * right.data[1].re - left.data[2].im * right.data[1].im
                        + left.data[3].re * right.data[3].re - left.data[3].im * right.data[3].im,
                left.data[2].re * right.data[1].im + left.data[2].im * right.data[1].re
                        + left.data[3].re * right.data[3].im + left.data[3].im * right.data[3].re);
    }

    private static void clearComplexMatrix(final ComplexMatrix matrix) {
        for (int i = 0; i < 4; ++i) {
            matrix.data[i].re = 0.0;
            matrix.data[i].im = 0.0;
        }
    }

    private static void addScaledComplexMatrix(final ComplexMatrix accumulator,
                                               final ComplexMatrix term,
                                               final double scaleRe,
                                               final double scaleIm) {
        for (int i = 0; i < 4; ++i) {
            final MutableComplex acc = accumulator.data[i];
            final MutableComplex value = term.data[i];
            final double re = scaleRe * value.re - scaleIm * value.im;
            final double im = scaleRe * value.im + scaleIm * value.re;
            acc.re += re;
            acc.im += im;
        }
    }

    private static void fillDividedDifferenceExpComplex(final double leftRe,
                                                        final double leftIm,
                                                        final double rightRe,
                                                        final double rightIm,
                                                        final MutableComplex out) {
        final double deltaRe = leftRe - rightRe;
        final double deltaIm = leftIm - rightIm;
        if (deltaRe * deltaRe + deltaIm * deltaIm < EPS * EPS) {
            fillExpComplex(leftRe, leftIm, out);
            return;
        }

        final double leftScale = Math.exp(leftRe);
        final double rightScale = Math.exp(rightRe);
        final double numeratorRe = leftScale * Math.cos(leftIm) - rightScale * Math.cos(rightIm);
        final double numeratorIm = leftScale * Math.sin(leftIm) - rightScale * Math.sin(rightIm);
        divideComplex(numeratorRe, numeratorIm, deltaRe, deltaIm, out);
    }

    private static void fillExpComplex(final double re,
                                       final double im,
                                       final MutableComplex out) {
        final double scale = Math.exp(re);
        out.re = scale * Math.cos(im);
        out.im = scale * Math.sin(im);
    }

    private static void setComplex(final MutableComplex out,
                                   final double re,
                                   final double im) {
        out.re = re;
        out.im = im;
    }

    private static final class PreparedCoefficients {
        private double alpha;
        private double beta;
        private double gamma;
        private double eta;
    }

    private static final class MutableSpectrum {
        private boolean distinct;
        private final MutableComplex lambda0 = new MutableComplex();
        private final MutableComplex lambda1 = new MutableComplex();
    }

    private static final class ComplexMatrix {
        private final MutableComplex[] data = new MutableComplex[]{
                new MutableComplex(),
                new MutableComplex(),
                new MutableComplex(),
                new MutableComplex()
        };
    }

    private static final class BuildWorkspace {
        private final double[] left = new double[4];
        private final double[] right = new double[4];
        private final double[] leftExp = new double[4];
        private final double[] rightExp = new double[4];
        private final double[] shifted = new double[4];
        private final double[] phi = new double[4];
        private final double[] basis = new double[4];
        private final double[] sylvester = new double[16];
        private final double[] rhs = new double[4];
        private final double[] solution = new double[4];
        private final double[] augmented = new double[20];
        private final double[] leftProduct = new double[4];
        private final double[] rightProduct = new double[4];
        private final double[] nilpotent = new double[4];
        private final double[] leftTerm = new double[4];
        private final double[] rightTerm = new double[4];
        private final double[] mixedTerm = new double[4];
        private final PreparedCoefficients prepared = new PreparedCoefficients();
        private final MutableSpectrum leftSpectrum = new MutableSpectrum();
        private final MutableSpectrum rightSpectrum = new MutableSpectrum();
        private final ComplexMatrix leftProjector0 = new ComplexMatrix();
        private final ComplexMatrix leftProjector1 = new ComplexMatrix();
        private final ComplexMatrix rightProjector0 = new ComplexMatrix();
        private final ComplexMatrix rightProjector1 = new ComplexMatrix();
        private final ComplexMatrix leftTimesBasis = new ComplexMatrix();
        private final ComplexMatrix complexTerm = new ComplexMatrix();
        private final ComplexMatrix complexAccumulator = new ComplexMatrix();
        private final MutableComplex complexScalar = new MutableComplex();
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
}
