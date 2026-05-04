package dr.evomodel.treedatalikelihood.continuous.backprop;

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;

/**
 * Pure block-space adjoint helper for Lyapunov contributions.
 */
public final class BlockDiagonalLyapunovAdjointHelper {

    private final int dim;
    private final BlockDiagonalLyapunovSolver solver;
    private final double[] blockDTransposeParams;
    private final Workspace ws;

    public BlockDiagonalLyapunovAdjointHelper(final int dim,
                                              final BlockDiagonalLyapunovSolver solver) {
        if (dim <= 0) {
            throw new IllegalArgumentException("dim must be positive");
        }
        this.dim = dim;
        this.solver = solver;
        this.blockDTransposeParams = new double[3 * dim - 2];
        this.ws = new Workspace(dim);
    }

    public void solveAdjointInDBasis(final DenseMatrix64F hBlock,
                                     final double[] blockDParams,
                                     final DenseMatrix64F yOut) {
        checkSquare(hBlock, "hBlock");
        checkSquare(yOut, "yOut");

        buildTransposeParams(blockDParams, blockDTransposeParams);
        copyNegated(hBlock, ws.negH);
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
        dst.set(src);
        CommonOps.scale(-1.0, dst);
    }

    private static void addDenseContribution(final DenseMatrix64F y,
                                             final DenseMatrix64F sigma,
                                             final DenseMatrix64F gradAcc,
                                             final DenseMatrix64F tmp) {
        MatrixOps.matMul(y.data, sigma.data, tmp.data, y.numRows);
        addEquals(gradAcc, tmp);
        multiplyTransposeLeft(y, sigma, tmp);
        addEquals(gradAcc, tmp);
    }

    private static void multiplyTransposeLeft(final DenseMatrix64F left,
                                              final DenseMatrix64F right,
                                              final DenseMatrix64F out) {
        final int dimension = left.numRows;
        final double[] leftData = left.data;
        final double[] rightData = right.data;
        final double[] outData = out.data;
        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += leftData[k * dimension + i] * rightData[k * dimension + j];
                }
                outData[rowOffset + j] = sum;
            }
        }
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
