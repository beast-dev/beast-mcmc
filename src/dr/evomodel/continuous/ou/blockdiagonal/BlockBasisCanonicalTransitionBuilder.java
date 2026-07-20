package dr.evomodel.continuous.ou.blockdiagonal;

import dr.evomodel.continuous.ou.math.OUPositiveDefiniteMatrixInversionOptions;
import dr.evomodel.continuous.ou.math.PositiveDefiniteMatrixOps;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalExpSolver;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import org.ejml.data.DenseMatrix64F;

/**
 * Rotated/block-basis OU work, Phase 2 (see
 * private-code/rotatedOU/rotated_ou_block_basis_plan.md).
 *
 * <p>Fills a {@link CanonicalGaussianTransition} entirely in the rotated
 * {@code Ytilde = T(Y - beta)} basis, for a single branch of length {@code dt}:
 *
 * <pre>
 *   Atilde_k     = e^{B dt}                                  (block-diagonal, compressed, O(p) storage)
 *   Omegatilde_k = Vtilde - Atilde_k Vtilde Atilde_k^T        (dense p x p, built in O(p^2) via
 *                                                               {@link BlockDiagonalMatrixOps}'s block-sparse
 *                                                               multiplies, exploiting Atilde_k's sparsity)
 * </pre>
 *
 * <p>Because the tree's working coordinate is centered at the stationary mean
 * ({@code Ytilde}'s own stationary mean is exactly zero by construction — see
 * the plan's Phase 2 section), every branch transition here is driftless:
 * {@code informationX = informationY = 0}. The four dense precision blocks
 * are then assembled from {@code Atilde_k} (block-sparse) and
 * {@code branchPrecision = Omegatilde_k^{-1}} (dense, via the same one
 * unavoidable {@code O(p^3)} positive-definite inversion any backend needs —
 * see the plan's section 1) using only block-sparse multiplies, so the only
 * per-branch {@code O(p^3)} cost left is that one inversion, not the two
 * dense rotations ({@code R*expD*R^{-1}} and {@code R*Omegatilde_k*R^T})
 * the existing natural-basis path pays every branch.</p>
 *
 * <p>{@code Vtilde} itself is supplied by a cached {@link BlockBasisGlobalSnapshot}
 * (Phase 1) and is not recomputed here.</p>
 */
public final class BlockBasisCanonicalTransitionBuilder {

    private final int dimension;
    private final BlockDiagonalExpSolver expSolver;
    private final double[] expD;
    private final DenseMatrix64F atildeVTilde;
    private final DenseMatrix64F omegaTilde;
    private final double[] omegaSymmetrizedFlat;
    private final double[] precisionFlat;
    private final double[] choleskyScratch;
    private final double[] lowerInverseScratch;
    private final double[] mFlat;
    private final double[] xxFlat;

    public BlockBasisCanonicalTransitionBuilder(final int dimension,
                                                final int[] blockStarts,
                                                final int[] blockSizes) {
        this.dimension = dimension;
        this.expSolver = new BlockDiagonalExpSolver(
                new BlockDiagonalExpSolver.BlockStructure(dimension, blockStarts, blockSizes));
        this.expD = new double[dimension + 2 * countTwoByTwoBlocks(blockSizes)];
        this.atildeVTilde = new DenseMatrix64F(dimension, dimension);
        this.omegaTilde = new DenseMatrix64F(dimension, dimension);
        this.omegaSymmetrizedFlat = new double[dimension * dimension];
        this.precisionFlat = new double[dimension * dimension];
        this.choleskyScratch = new double[dimension * dimension];
        this.lowerInverseScratch = new double[dimension * dimension];
        this.mFlat = new double[dimension * dimension];
        this.xxFlat = new double[dimension * dimension];
    }

    public void fillCanonicalTransition(final BlockBasisGlobalSnapshot snapshot,
                                        final double dt,
                                        final CanonicalGaussianTransition out) {
        snapshot.ensureCurrent();
        final int[] blockStarts = snapshot.getBlockStarts();
        final int[] blockSizes = snapshot.getBlockSizes();
        final DenseMatrix64F vTilde = snapshot.getStationaryCovarianceDBasis();

        expSolver.computeCompressed(snapshot.getBlockDiagonalCompressed(), dt, expD);

        // Omegatilde = Vtilde - Atilde * Vtilde * Atilde^T, both multiplies O(p^2)
        // thanks to Atilde's block-diagonal sparsity.
        BlockDiagonalMatrixOps.multiplyBlockDiagonalLeft(
                expD, vTilde.data, atildeVTilde.data, dimension, blockStarts, blockSizes);
        BlockDiagonalMatrixOps.multiplyRightBlockDiagonalTranspose(
                atildeVTilde.data, expD, omegaTilde.data, dimension, blockStarts, blockSizes);
        BlockDiagonalDenseMatrixOps.subtract(vTilde, omegaTilde, omegaTilde);
        symmetrize(omegaTilde);

        final double logDet = assemblePrecisionBlocksFromInverse(blockStarts, blockSizes);

        // Atilde^T * M = Atilde^T * branchPrecision * Atilde (O(p^2), Atilde sparse).
        BlockDiagonalMatrixOps.multiplyBlockDiagonalLeftTranspose(
                expD, mFlat, xxFlat, dimension, blockStarts, blockSizes);

        final int d2 = dimension * dimension;
        System.arraycopy(precisionFlat, 0, out.precisionYY, 0, d2);
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                out.precisionYX[i * dimension + j] = -mFlat[i * dimension + j];
                out.precisionXY[j * dimension + i] = -mFlat[i * dimension + j];
            }
        }
        System.arraycopy(xxFlat, 0, out.precisionXX, 0, d2);
        symmetrizeFlat(out.precisionXX, dimension);

        java.util.Arrays.fill(out.informationX, 0, dimension, 0.0);
        java.util.Arrays.fill(out.informationY, 0, dimension, 0.0);

        out.logNormalizer = 0.5 * (dimension * Math.log(2.0 * Math.PI) + logDet);
    }

    private double assemblePrecisionBlocksFromInverse(final int[] blockStarts,
                                                      final int[] blockSizes) {
        // The one unavoidable O(p^3) step: factor and invert the dense branch covariance.
        final double logDet = PositiveDefiniteMatrixOps.copyAndInvertPositiveDefiniteFlat(
                omegaTilde, omegaSymmetrizedFlat, precisionFlat,
                choleskyScratch, lowerInverseScratch,
                OUPositiveDefiniteMatrixInversionOptions.OU_TREE);

        // M = branchPrecision * Atilde (O(p^2), Atilde sparse).
        BlockDiagonalMatrixOps.multiplyRightBlockDiagonal(
                precisionFlat, expD, dimension, blockStarts, blockSizes, mFlat);
        return logDet;
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

    private static int countTwoByTwoBlocks(final int[] blockSizes) {
        int count = 0;
        for (final int size : blockSizes) {
            if (size == 2) {
                count++;
            }
        }
        return count;
    }

    private static void symmetrizeFlat(final double[] data, final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = i + 1; j < dimension; ++j) {
                final double value = 0.5 * (data[i * dimension + j] + data[j * dimension + i]);
                data[i * dimension + j] = value;
                data[j * dimension + i] = value;
            }
        }
    }
}
