package dr.evomodel.treedatalikelihood.continuous.backprop;

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.Arrays;

/**
 * Pure block-space adjoint helper for Lyapunov contributions.
 */
public final class BlockDiagonalLyapunovAdjointHelper {

    private static final String USE_EQUAL_DIAGONAL_PLAN_PROPERTY =
            "beast.experimental.useEqualDiagonalLyapunovPlan";
    private static final boolean USE_EQUAL_DIAGONAL_PLAN =
            Boolean.parseBoolean(System.getProperty(USE_EQUAL_DIAGONAL_PLAN_PROPERTY, "true"));
    private static final double EQUAL_DIAGONAL_TOLERANCE = 1.0e-12;

    private final int dim;
    private final BlockDiagonalLyapunovSolver solver;
    private final double[] blockDTransposeParams;
    private final boolean equalDiagonalStructureCandidate;
    private final BlockDiagonalLyapunovEqualDiagonalPlan.Plan equalDiagonalPlan;
    private final int[] equalDiagonalBlock2Starts;
    private final double[] cachedEqualDiagonalPlanParams;
    private int cachedEqualDiagonalPlanParameterHash;
    private boolean cachedEqualDiagonalPlanValid;
    private final Workspace ws;

    public BlockDiagonalLyapunovAdjointHelper(final int dim,
                                              final BlockDiagonalLyapunovSolver solver) {
        this(dim, solver, null);
    }

    public BlockDiagonalLyapunovAdjointHelper(final int dim,
                                              final BlockDiagonalLyapunovSolver solver,
                                              final BlockDiagonalExpSolver.BlockStructure structure) {
        if (dim <= 0) {
            throw new IllegalArgumentException("dim must be positive");
        }
        this.dim = dim;
        this.solver = solver;
        this.blockDTransposeParams = new double[3 * dim - 2];
        this.equalDiagonalStructureCandidate = USE_EQUAL_DIAGONAL_PLAN
                && structure != null
                && isLeadingScalarPlusTwoByTwoStructure(structure);
        this.equalDiagonalPlan = equalDiagonalStructureCandidate
                ? BlockDiagonalLyapunovEqualDiagonalPlan.createPlan(structure)
                : null;
        this.equalDiagonalBlock2Starts = equalDiagonalStructureCandidate
                ? collectTwoByTwoBlockStarts(structure)
                : null;
        this.cachedEqualDiagonalPlanParams = equalDiagonalStructureCandidate ? new double[3 * dim - 2] : null;
        this.cachedEqualDiagonalPlanParameterHash = 0;
        this.cachedEqualDiagonalPlanValid = false;
        this.ws = new Workspace(dim);
    }

    public void solveAdjointInDBasis(final DenseMatrix64F hBlock,
                                     final double[] blockDParams,
                                     final DenseMatrix64F yOut) {
        checkSquare(hBlock, "hBlock");
        checkSquare(yOut, "yOut");

        copyNegated(hBlock, ws.negH);
        if (trySolveWithEqualDiagonalPlan(blockDParams, ws.negH, yOut)) {
            return;
        }
        buildTransposeParams(blockDParams, blockDTransposeParams);
        solver.solve(blockDTransposeParams, ws.negH, yOut);
    }

    public void accumulateLyapunovContributionInDBasis(final DenseMatrix64F sigmaBlock,
                                                       final DenseMatrix64F hBlock,
                                                       final double[] blockDParams,
                                                       final DenseMatrix64F gradDAcc) {
        checkSquare(sigmaBlock, "sigmaBlock");
        checkSquare(hBlock, "hBlock");
        checkSquare(gradDAcc, "gradDAcc");

        solveAdjointInDBasis(hBlock, blockDParams, ws.y);
        addDenseContribution(ws.y, sigmaBlock, gradDAcc, ws.tmp);
    }

    private void buildTransposeParams(final double[] src, final double[] dst) {
        final int offLen = dim - 1;
        System.arraycopy(src, 0, dst, 0, dim);
        System.arraycopy(src, dim + offLen, dst, dim, offLen);
        System.arraycopy(src, dim, dst, dim + offLen, offLen);
    }

    private static void copyNegated(final DenseMatrix64F src, final DenseMatrix64F dst) {
        final int length = src.numRows * src.numCols;
        final double[] srcData = src.data;
        final double[] dstData = dst.data;
        for (int i = 0; i < length; ++i) {
            dstData[i] = -srcData[i];
        }
    }

    private boolean trySolveWithEqualDiagonalPlan(final double[] blockDParams,
                                                 final DenseMatrix64F negH,
                                                 final DenseMatrix64F yOut) {
        // The cached plan mirrors upper-block solves, so only use it when the
        // adjoint right-hand side is symmetric. Generic block solves handle the
        // full non-symmetric case.
        if (!equalDiagonalStructureCandidate
                || !isSymmetric(negH)
                || !hasEqualDiagonalParameters(blockDParams)) {
            return false;
        }

        final int parameterHash = Arrays.hashCode(blockDParams);
        final boolean parametersChanged = !cachedEqualDiagonalPlanValid
                || parameterHash != cachedEqualDiagonalPlanParameterHash
                || !Arrays.equals(blockDParams, cachedEqualDiagonalPlanParams);

        if (parametersChanged) {
            System.arraycopy(blockDParams, 0, cachedEqualDiagonalPlanParams, 0, blockDParams.length);
            cachedEqualDiagonalPlanParameterHash = parameterHash;
            cachedEqualDiagonalPlanValid = equalDiagonalPlan.updateTransposedParameters(blockDParams);
            if (!cachedEqualDiagonalPlanValid) {
                return false;
            }
        }

        return equalDiagonalPlan.tryApplySymmetric(negH, yOut);
    }

    private boolean isSymmetric(final DenseMatrix64F matrix) {
        final double[] data = matrix.data;
        for (int i = 0; i < dim; ++i) {
            final int rowOffset = i * dim;
            for (int j = i + 1; j < dim; ++j) {
                if (Math.abs(data[rowOffset + j] - data[j * dim + i]) > EQUAL_DIAGONAL_TOLERANCE) {
                    return false;
                }
            }
        }
        return true;
    }

    private boolean hasEqualDiagonalParameters(final double[] blockDParams) {
        final int[] starts = equalDiagonalBlock2Starts;
        for (int k = 0; k < starts.length; ++k) {
            final int start = starts[k];
            if (Math.abs(blockDParams[start] - blockDParams[start + 1]) > EQUAL_DIAGONAL_TOLERANCE) {
                return false;
            }
        }
        return true;
    }

    private static int[] collectTwoByTwoBlockStarts(final BlockDiagonalExpSolver.BlockStructure structure) {
        int count = 0;
        for (int b = 0; b < structure.getNumBlocks(); ++b) {
            if (structure.getBlockSize(b) == 2) {
                count++;
            }
        }
        final int[] starts = new int[count];
        for (int b = 0, k = 0; b < structure.getNumBlocks(); ++b) {
            if (structure.getBlockSize(b) == 2) {
                starts[k++] = structure.getBlockStart(b);
            }
        }
        return starts;
    }

    private static boolean isLeadingScalarPlusTwoByTwoStructure(
            final BlockDiagonalExpSolver.BlockStructure structure) {
        final int blockCount = structure.getNumBlocks();
        if (blockCount == 0) {
            return false;
        }
        int block = 0;
        if (structure.getBlockSize(0) == 1) {
            block = 1;
        }
        for (; block < blockCount; ++block) {
            if (structure.getBlockSize(block) != 2) {
                return false;
            }
        }
        return true;
    }

    private static void addDenseContribution(final DenseMatrix64F y,
                                             final DenseMatrix64F sigma,
                                             final DenseMatrix64F gradAcc,
                                             final DenseMatrix64F tmp) {
        CommonOps.mult(y, sigma, tmp);
        addEquals(gradAcc, tmp);
        multiplyTransposeLeft(y, sigma, tmp);
        addEquals(gradAcc, tmp);
    }

    private static void multiplyTransposeLeft(final DenseMatrix64F left,
                                              final DenseMatrix64F right,
                                              final DenseMatrix64F out) {
        CommonOps.multTransA(left, right, out);
    }

    private static void addEquals(final DenseMatrix64F target,
                                  final DenseMatrix64F increment) {
        final int length = target.numRows * target.numCols;
        final double[] targetData = target.data;
        final double[] incrementData = increment.data;
        for (int i = 0; i < length; ++i) {
            targetData[i] += incrementData[i];
        }
    }

    private void checkSquare(final DenseMatrix64F m, final String name) {
        if (m.numRows != dim || m.numCols != dim) {
            throw new IllegalArgumentException(name + " must be " + dim + " x " + dim);
        }
    }

    private static final class Workspace {
        final DenseMatrix64F negH;
        final DenseMatrix64F y;
        final DenseMatrix64F tmp;

        private Workspace(final int dim) {
            this.negH = new DenseMatrix64F(dim, dim);
            this.y = new DenseMatrix64F(dim, dim);
            this.tmp = new DenseMatrix64F(dim, dim);
        }
    }
}
