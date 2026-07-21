package dr.evomodel.continuous.ou.blockdiagonal;

import dr.evomodel.continuous.ou.blockdiagonal.math.BlockDiagonalExpSolver;
import dr.evomodel.continuous.ou.blockdiagonal.math.BlockDiagonalFrechetHelper;
import dr.evomodel.continuous.ou.blockdiagonal.math.BlockDiagonalLyapunovAdjointHelper;
import dr.evomodel.continuous.ou.blockdiagonal.math.BlockDiagonalLyapunovSolver;
import dr.evomodel.continuous.canonical.math.MatrixOps;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalLocalTransitionAdjoints;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * Converts per-branch {@code Ytilde}-basis transition adjoints
 * ({@code dLogL_dF} where {@code F = Atilde_k}, {@code dLogL_dOmega} where
 * {@code Omega = Omegatilde_k}) into gradients with respect to the
 * compressed block-diagonal {@code B} parameters, the rotation {@code R},
 * and the diffusion matrix {@code Q}.
 *
 * <p>This mirrors {@link BlockDiagonalSelectionMatrixParameterization}'s
 * sandwich-adjoint code, but simpler: because the new backend builds
 * {@code Atilde_k}/{@code Omegatilde_k} directly with no per-branch
 * rotation, {@code dLogL_dF}/{@code dLogL_dOmega} are already expressed
 * w.r.t. the D-basis quantities themselves, so none of the dense backend's
 * {@code R (...) R^{-1}}/{@code R (...) R^T} sandwiches are needed here.
 * The one thing that <em>is</em> shared/global rather than per-branch is
 * {@code Vtilde} itself (built once by {@link BlockBasisGlobalSnapshot}, not
 * per branch), so each branch's "direct" contribution to {@code dLogL_dVtilde}
 * is accumulated into a caller-supplied running total and only converted to
 * {@code R}/{@code Q} gradients once, via {@link #finalizeGlobalPullback}.
 *
 * <p>Stationary-mean ({@code beta}) gradients are out of scope for this
 * class.
 *
 * <p>Status: not yet wired into the production tree-likelihood traversal;
 * see {@link BlockBasisGlobalSnapshot} for the current status and the
 * wiring plan.</p>
 */
public final class BlockBasisCanonicalAdjointBuilder {

    private final BlockDiagonalStructure structure;
    private final int dimension;
    private final int[] blockStarts;
    private final int[] blockSizes;
    private final BlockDiagonalExpSolver expSolver;
    private final BlockDiagonalFrechetHelper frechetHelper;
    private final BlockDiagonalLyapunovAdjointHelper lyapunovAdjointHelper;

    private final double[] expD;
    private final DenseMatrix64F dLogL_dF_mat;
    private final DenseMatrix64F dLogL_dOmega_mat;
    private final DenseMatrix64F gradD;
    private final DenseMatrix64F temp1;
    private final DenseMatrix64F temp2;
    private final DenseMatrix64F gECov;
    private final DenseMatrix64F yAdjoint;
    private final DenseMatrix64F negatedGS;

    public BlockBasisCanonicalAdjointBuilder(final BlockDiagonalStructure structure) {
        this.structure = structure;
        this.dimension = structure.getDimension();
        this.blockStarts = structure.getBlockStarts();
        this.blockSizes = structure.getBlockSizes();
        this.expSolver = new BlockDiagonalExpSolver(structure.getMathStructure());
        this.frechetHelper = new BlockDiagonalFrechetHelper(structure.getMathStructure());
        final BlockDiagonalLyapunovSolver solver =
                new BlockDiagonalLyapunovSolver(dimension, blockStarts, blockSizes);
        this.lyapunovAdjointHelper =
                new BlockDiagonalLyapunovAdjointHelper(dimension, solver, structure.getMathStructure());

        this.expD = new double[structure.getCompressedDimension()];
        this.dLogL_dF_mat = new DenseMatrix64F(dimension, dimension);
        this.dLogL_dOmega_mat = new DenseMatrix64F(dimension, dimension);
        this.gradD = new DenseMatrix64F(dimension, dimension);
        this.temp1 = new DenseMatrix64F(dimension, dimension);
        this.temp2 = new DenseMatrix64F(dimension, dimension);
        this.gECov = new DenseMatrix64F(dimension, dimension);
        this.yAdjoint = new DenseMatrix64F(dimension, dimension);
        this.negatedGS = new DenseMatrix64F(dimension, dimension);
    }

    /**
     * Accumulates one branch's contribution into {@code compressedDAccumulator}
     * (active-compressed format, size {@code dimension + 2*num2x2Blocks}) and
     * into the global {@code gSAccumulator} (dense, {@code dLogL_dVtilde}'s
     * "direct" part, summed across all branches before being finalized).
     */
    public void accumulateBranchPullback(final BlockBasisGlobalSnapshot snapshot,
                                         final double dt,
                                         final CanonicalLocalTransitionAdjoints adjoints,
                                         final double[] compressedDAccumulator,
                                         final DenseMatrix64F gSAccumulator) {
        snapshot.ensureCurrent();
        structure.checkMatches(snapshot.getStructure());
        final double[] blockDParams = snapshot.getBlockDiagonalCompressed();
        final DenseMatrix64F vTilde = snapshot.getStationaryCovarianceDBasis();
        expSolver.computeCompressed(blockDParams, dt, expD);

        // F = Atilde_k adjoint -> compressed-B gradient. No R-sandwich: F is Atilde_k directly.
        System.arraycopy(adjoints.dLogL_dF, 0, dLogL_dF_mat.data, 0, dimension * dimension);
        frechetHelper.frechetAdjointExpInDBasis(blockDParams, dLogL_dF_mat, dt, gradD);
        transposeInPlace(gradD);
        accumulateCompressedGradient(gradD, compressedDAccumulator);

        // Omega = Vtilde - Atilde_k Vtilde Atilde_k^T adjoint, decomposed into:
        // (a) the "direct" Vtilde-dependence part, folded into the global gS total, and
        // (b) the Atilde_k-dependence part (Frechet, per branch since Atilde_k is dt-specific).
        System.arraycopy(adjoints.dLogL_dOmega, 0, dLogL_dOmega_mat.data, 0, dimension * dimension);
        symmetrize(dLogL_dOmega_mat);

        BlockDiagonalMatrixOps.multiplyBlockDiagonalLeftTranspose(
                expD, dLogL_dOmega_mat.data, temp1.data, dimension, blockStarts, blockSizes);
        BlockDiagonalMatrixOps.multiplyRightBlockDiagonal(
                temp1.data, expD, dimension, blockStarts, blockSizes, temp2.data);
        BlockDiagonalDenseMatrixOps.subtract(dLogL_dOmega_mat, temp2, temp2);
        addDenseMatrixToFlatArray(temp2, gSAccumulator);

        BlockDiagonalMatrixOps.multiplyRightBlockDiagonal(
                dLogL_dOmega_mat.data, expD, dimension, blockStarts, blockSizes, temp1.data);
        CommonOps.mult(temp1, vTilde, gECov);
        CommonOps.scale(-2.0, gECov);
        frechetHelper.frechetAdjointExpInDBasis(blockDParams, gECov, dt, gradD);
        transposeInPlace(gradD);
        accumulateCompressedGradient(gradD, compressedDAccumulator);
    }

    /**
     * Folds the root prior's direct dependence on {@code Vtilde} into the same
     * global total as {@link #accumulateBranchPullback}'s "direct" part (the
     * root prior covariance is {@code Vtilde} itself, with no {@code Atilde_k}
     * sandwich at all, so there is no Frechet/B-parameter contribution here).
     */
    public void accumulateRootPriorPullback(final double[] dLogL_dPriorCovFlat,
                                            final DenseMatrix64F gSAccumulator) {
        final int d2 = dimension * dimension;
        final double[] data = gSAccumulator.data;
        for (int i = 0; i < d2; ++i) {
            data[i] += dLogL_dPriorCovFlat[i];
        }
    }

    /**
     * One-time global step: propagates the accumulated {@code Vtilde}-adjoint
     * total to the compressed-B gradient and to {@code R}/{@code Q} via a
     * single Lyapunov-adjoint solve (the "accumulated globally rather than
     * per branch" step the plan calls for).
     */
    public void finalizeGlobalPullback(final BlockBasisGlobalSnapshot snapshot,
                                       final DenseMatrix64F gSAccumulator,
                                       final double[] compressedDAccumulator,
                                       final DenseMatrix64F rGradAccumulator,
                                       final double[] qGradFlatAccumulator) {
        snapshot.ensureCurrent();
        structure.checkMatches(snapshot.getStructure());
        final double[] blockDParams = snapshot.getBlockDiagonalCompressed();
        symmetrize(gSAccumulator);

        CommonOps.fill(gradD, 0.0);
        lyapunovAdjointHelper.accumulateLyapunovContributionInDBasis(
                snapshot.getStationaryCovarianceDBasis(), gSAccumulator, blockDParams, gradD);
        accumulateCompressedGradient(gradD, compressedDAccumulator);

        lyapunovAdjointHelper.solveAdjointInDBasis(gSAccumulator, blockDParams, yAdjoint);

        final DenseMatrix64F qTilde = snapshot.getQDBasis();
        CommonOps.mult(yAdjoint, qTilde, temp1);
        CommonOps.multTransA(yAdjoint, qTilde, temp2);
        CommonOps.add(temp1, temp2, gECov);
        CommonOps.multTransA(snapshot.getRotationInverse(), gECov, temp1);
        addDenseMatrixToFlatArray(temp1, rGradAccumulator);

        negatedGS.set(gSAccumulator);
        CommonOps.scale(-1.0, negatedGS);
        lyapunovAdjointHelper.solveAdjointInDBasis(negatedGS, blockDParams, yAdjoint);
        MatrixOps.symmetricSandwichTransposeLeft(
                snapshot.getRotationInverse().data, yAdjoint.data, temp1.data, temp2.data, dimension);
        for (int i = 0; i < dimension * dimension; ++i) {
            qGradFlatAccumulator[i] += temp1.data[i];
        }
    }

    /**
     * Explicit Jacobian-term gradient contribution: {@code Rbar += -N_obs * R^{-T}},
     * where {@code R^{-T} = T^T} ({@code T = R^{-1}}).
     */
    public void addJacobianTermGradient(final BlockBasisGlobalSnapshot snapshot,
                                        final int numObservedTips,
                                        final DenseMatrix64F rGradAccumulator) {
        structure.checkMatches(snapshot.getStructure());
        final DenseMatrix64F t = snapshot.getRotationInverse();
        final double[] data = rGradAccumulator.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                data[i * dimension + j] += -numObservedTips * t.get(j, i);
            }
        }
    }

    /**
     * One tip's {@code R}-gradient contribution from {@code ytilde = T(y - beta)}'s own
     * dependence on {@code R} (via {@code T = R^{-1}}), holding {@code y}/{@code beta} fixed:
     * {@code Rbar += -T^T * dLogL_dYTilde * ytilde^T}.
     */
    public void addTipTransformRGradient(final BlockBasisGlobalSnapshot snapshot,
                                        final double[] yTilde,
                                        final double[] dLogL_dYTilde,
                                        final DenseMatrix64F rGradAccumulator) {
        structure.checkMatches(snapshot.getStructure());
        final DenseMatrix64F t = snapshot.getRotationInverse();
        final double[] data = rGradAccumulator.data;
        for (int i = 0; i < dimension; ++i) {
            double tTg = 0.0;
            for (int k = 0; k < dimension; ++k) {
                tTg += t.get(k, i) * dLogL_dYTilde[k];
            }
            for (int j = 0; j < dimension; ++j) {
                data[i * dimension + j] += -tTg * yTilde[j];
            }
        }
    }

    private void accumulateCompressedGradient(final DenseMatrix64F denseGradient,
                                              final double[] compressedAccumulator) {
        final double[] data = denseGradient.data;
        final int upperBase = dimension;
        final int lowerBase = dimension + structure.getNumTwoByTwoBlocks();
        int blockIndex = 0;

        for (int i = 0; i < dimension; ++i) {
            compressedAccumulator[i] += data[i * dimension + i];
        }
        for (int b = 0; b < blockStarts.length; ++b) {
            if (blockSizes[b] != 2) {
                continue;
            }
            final int start = blockStarts[b];
            compressedAccumulator[upperBase + blockIndex] += data[start * dimension + start + 1];
            compressedAccumulator[lowerBase + blockIndex] += data[(start + 1) * dimension + start];
            blockIndex++;
        }
    }

    private static void addDenseMatrixToFlatArray(final DenseMatrix64F src, final DenseMatrix64F dest) {
        addDenseMatrixToFlatArray(src, dest.data);
    }

    private static void addDenseMatrixToFlatArray(final DenseMatrix64F src, final double[] dest) {
        final int length = src.numRows * src.numCols;
        final double[] data = src.data;
        for (int i = 0; i < length; ++i) {
            dest[i] += data[i];
        }
    }

    private static void symmetrize(final DenseMatrix64F matrix) {
        final int dim = matrix.numRows;
        final double[] data = matrix.data;
        for (int i = 0; i < dim; ++i) {
            for (int j = i + 1; j < dim; ++j) {
                final double value = 0.5 * (data[i * dim + j] + data[j * dim + i]);
                data[i * dim + j] = value;
                data[j * dim + i] = value;
            }
        }
    }

    private static void transposeInPlace(final DenseMatrix64F matrix) {
        final int dim = matrix.numRows;
        final double[] data = matrix.data;
        for (int i = 0; i < dim; ++i) {
            for (int j = i + 1; j < dim; ++j) {
                final int ij = i * dim + j;
                final int ji = j * dim + i;
                final double tmp = data[ij];
                data[ij] = data[ji];
                data[ji] = tmp;
            }
        }
    }

}
