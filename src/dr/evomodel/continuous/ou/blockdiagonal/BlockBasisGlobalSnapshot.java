package dr.evomodel.continuous.ou.blockdiagonal;

import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalLyapunovSolver;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.BlockDiagonalDecomposition;
import dr.inference.model.MatrixParameterInterface;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.Arrays;

/**
 * Rotated/block-basis OU work (see
 * private-code/rotatedOU/rotated_ou_block_basis_plan.md, Phase 1):
 * one-time, per-model-change cache of the quantities a block-basis
 * evaluation path needs but that do not depend on branch length:
 *
 * <pre>
 *   T            = R^{-1}            (change of basis to Ytilde = T(Y - beta))
 *   Vtilde       = T V T^T           (stationary covariance in the D-basis,
 *                                      via the block Lyapunov solve, not a
 *                                      dense p x p sandwich)
 *   log|det R|   = Jacobian term for the change of variables
 * </pre>
 *
 * <p>This intentionally re-derives R/R^{-1}/Vtilde independently of
 * {@link BlockDiagonalSelectionMatrixParameterization} (which recomputes the
 * D-basis stationary covariance inside every per-branch
 * {@code fillTransitionCovarianceMatrix} call) so a single instance here can
 * be cached once per likelihood evaluation and reused across all branches,
 * instead of being paid {@code O(p^3)} per branch.</p>
 */
public final class BlockBasisGlobalSnapshot {

    private final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter;
    private final MatrixParameterInterface diffusionMatrix;
    private final int dimension;
    private final BlockDiagonalDecomposition decomposition;
    private final BlockDiagonalLyapunovSolver lyapunovSolver;

    private final DenseMatrix64F rMatrix;
    private final DenseMatrix64F rinvMatrix;
    private final DenseMatrix64F qMatrix;
    private final DenseMatrix64F rinvQ;
    private final DenseMatrix64F qDBasis;
    private final DenseMatrix64F stationaryCovDBasis;

    private double logAbsDetR;

    private boolean valid;
    private final double[] cachedR;
    private final double[] cachedRinv;
    private final double[] cachedBlockD;
    private final double[] cachedQ;

    public BlockBasisGlobalSnapshot(final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                    final MatrixParameterInterface diffusionMatrix) {
        this.blockParameter = blockParameter;
        this.diffusionMatrix = diffusionMatrix;
        this.dimension = blockParameter.getRowDimension();
        this.decomposition = blockParameter.createBlockDiagonalDecomposition();
        this.lyapunovSolver = new BlockDiagonalLyapunovSolver(
                dimension, blockParameter.getBlockStarts(), blockParameter.getBlockSizes());

        this.rMatrix = new DenseMatrix64F(dimension, dimension);
        this.rinvMatrix = new DenseMatrix64F(dimension, dimension);
        this.qMatrix = new DenseMatrix64F(dimension, dimension);
        this.rinvQ = new DenseMatrix64F(dimension, dimension);
        this.qDBasis = new DenseMatrix64F(dimension, dimension);
        this.stationaryCovDBasis = new DenseMatrix64F(dimension, dimension);

        this.cachedR = new double[dimension * dimension];
        this.cachedRinv = new double[dimension * dimension];
        this.cachedBlockD = new double[blockParameter.getTridiagonalDDimension()];
        this.cachedQ = new double[dimension * dimension];
        Arrays.fill(cachedR, Double.NaN);
        this.valid = false;
    }

    public int getDimension() {
        return dimension;
    }

    /** Refreshes {@code R}, {@code T = R^{-1}}, {@code Vtilde}, and {@code log|det R|} if anything changed. */
    public void ensureCurrent() {
        blockParameter.fillBlockDiagonalDecomposition(decomposition);
        final double[] rData = decomposition.getR();
        final double[] rinvData = decomposition.getRInverse();
        final double[] blockD = decomposition.getBlockDiagonal();

        if (valid
                && Arrays.equals(rData, cachedR)
                && Arrays.equals(rinvData, cachedRinv)
                && Arrays.equals(blockD, cachedBlockD)
                && diffusionMatches()) {
            return;
        }

        System.arraycopy(rData, 0, rMatrix.data, 0, rData.length);
        System.arraycopy(rinvData, 0, rinvMatrix.data, 0, rinvData.length);
        fillDenseMatrix(diffusionMatrix, qMatrix);

        // qDBasis = T Q T^T, with T = R^{-1}.
        CommonOps.mult(rinvMatrix, qMatrix, rinvQ);
        MatrixOps.matMulTransposedRight(rinvQ.data, rinvMatrix.data, qDBasis.data, dimension);
        symmetrize(qDBasis);

        // Vtilde solves B*Vtilde + Vtilde*B^T + qDBasis = 0 in the D-basis directly,
        // never materializing a dense p x p T*V*T^T sandwich.
        lyapunovSolver.solve(blockD, qDBasis, stationaryCovDBasis);

        logAbsDetR = Math.log(Math.abs(CommonOps.det(rMatrix)));

        System.arraycopy(rData, 0, cachedR, 0, rData.length);
        System.arraycopy(rinvData, 0, cachedRinv, 0, rinvData.length);
        System.arraycopy(blockD, 0, cachedBlockD, 0, blockD.length);
        cacheDiffusion();
        valid = true;
    }

    /** {@code R}, the rotation/change-of-basis matrix in {@code S = R B R^{-1}}. Read-only. */
    public DenseMatrix64F getRotation() {
        return rMatrix;
    }

    /** {@code T = R^{-1}}, mapping natural coordinates to the rotated basis {@code Ytilde = T(Y - beta)}. Read-only. */
    public DenseMatrix64F getRotationInverse() {
        return rinvMatrix;
    }

    /** {@code Vtilde = T V T^T}, the stationary covariance expressed in the D-basis. Read-only. */
    public DenseMatrix64F getStationaryCovarianceDBasis() {
        return stationaryCovDBasis;
    }

    public double getLogAbsDetR() {
        return logAbsDetR;
    }

    /** {@code Qtilde = T Q T^T}, the diffusion matrix expressed in the D-basis. Read-only. */
    public DenseMatrix64F getQDBasis() {
        return qDBasis;
    }

    /** Compressed tridiagonal-form {@code B} (diag/upper/lower), as used by {@link BlockDiagonalExpSolver}. Read-only. */
    public double[] getBlockDiagonalCompressed() {
        return decomposition.getBlockDiagonal();
    }

    public int[] getBlockStarts() {
        return decomposition.getBlockStarts();
    }

    public int[] getBlockSizes() {
        return decomposition.getBlockSizes();
    }

    private boolean diffusionMatches() {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                if (diffusionMatrix.getParameterValue(i, j) != cachedQ[i * dimension + j]) {
                    return false;
                }
            }
        }
        return true;
    }

    private void cacheDiffusion() {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                cachedQ[i * dimension + j] = diffusionMatrix.getParameterValue(i, j);
            }
        }
    }

    private static void fillDenseMatrix(final MatrixParameterInterface parameter, final DenseMatrix64F out) {
        final int dimension = out.numRows;
        final double[] data = out.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                data[i * dimension + j] = parameter.getParameterValue(i, j);
            }
        }
    }

    private static void symmetrize(final DenseMatrix64F matrix) {
        final int dimension = matrix.numRows;
        final double[] data = matrix.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = i + 1; j < dimension; ++j) {
                final double value = 0.5 * (data[i * dimension + j] + data[j * dimension + i]);
                data[i * dimension + j] = value;
                data[j * dimension + i] = value;
            }
        }
    }
}
