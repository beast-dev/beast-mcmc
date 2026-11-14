package dr.evomodel.continuous.ou.blockdiagonal;

import dr.evomodel.continuous.ou.math.OUPositiveDefiniteMatrixInversionOptions;
import dr.evomodel.continuous.ou.math.PositiveDefiniteMatrixOps;
import dr.evomodel.continuous.ou.blockdiagonal.math.BlockDiagonalExpSolver;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianTransition;
import org.ejml.data.DenseMatrix64F;

import java.util.Arrays;

/**
 * Fills a {@link CanonicalGaussianTransition} entirely in the rotated
 * {@code Ytilde = T(Y - beta)} basis, for a single branch of length {@code dt}:
 *
 * <pre>
 *   Atilde_k     = e^{-B dt}                                 (block-diagonal, compressed, O(p) storage)
 *   Omegatilde_k = Vtilde - Atilde_k Vtilde Atilde_k^T        (dense p x p, built in O(p^2) via
 *                                                               {@link BlockDiagonalMatrixOps}'s block-sparse
 *                                                               multiplies, exploiting Atilde_k's sparsity)
 * </pre>
 *
 * <p>Because the tree's working coordinate is centered at the stationary mean
 * ({@code Ytilde}'s own stationary mean is exactly zero by construction),
 * every branch transition here is driftless:
 * {@code informationX = informationY = 0}. The four dense precision blocks
 * are then assembled from {@code Atilde_k} (block-sparse) and
 * {@code branchPrecision = Omegatilde_k^{-1}} (dense, via the same one
 * unavoidable {@code O(p^3)} positive-definite inversion any backend needs)
 * using only block-sparse multiplies, so the only per-branch {@code O(p^3)}
 * cost left is that one inversion, not the two dense rotations
 * ({@code R*expD*R^{-1}} and {@code R*Omegatilde_k*R^T}) the existing
 * natural-basis path pays every branch.</p>
 *
 * <p>{@code Vtilde} itself is supplied by a cached {@link BlockBasisGlobalSnapshot}
 * and is not recomputed here.</p>
 */
public final class BlockBasisCanonicalTransitionBuilder {

    private final BlockDiagonalStructure structure;
    private final int dimension;
    private final BlockDiagonalExpSolver expSolver;
    private final double[] expD;
    private final double[] dExpDdt;
    private final DenseMatrix64F atildeVTilde;
    private final DenseMatrix64F omegaTilde;
    private final double[] omegaSymmetrizedFlat;
    private final double[] precisionFlat;
    private final double[] choleskyScratch;
    private final double[] lowerInverseScratch;
    private final double[] mFlat;
    private final double[] xxFlat;

    public BlockBasisCanonicalTransitionBuilder(final BlockDiagonalStructure structure) {
        this.structure = structure;
        final int dimension = structure.getDimension();
        this.dimension = dimension;
        this.expSolver = new BlockDiagonalExpSolver(structure.getMathStructure());
        this.expD = new double[structure.getCompressedDimension()];
        this.dExpDdt = new double[structure.getCompressedDimension()];
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
        computeTransitionMoments(snapshot, dt);
        fillCanonicalTransitionFromCurrentMoments(snapshot, out);
    }

    public void fillCanonicalTransitionAndMoments(final BlockBasisGlobalSnapshot snapshot,
                                                  final double dt,
                                                  final CanonicalGaussianTransition out,
                                                  final double[] transitionMatrixOut,
                                                  final double[] covarianceOut) {
        computeTransitionMoments(snapshot, dt);
        fillTransitionMatrixFromCompressed(
                expD, transitionMatrixOut, dimension, snapshot.getBlockStarts(), snapshot.getBlockSizes());
        System.arraycopy(omegaTilde.data, 0, covarianceOut, 0, dimension * dimension);
        fillCanonicalTransitionFromCurrentMoments(snapshot, out);
    }

    private void fillCanonicalTransitionFromCurrentMoments(final BlockBasisGlobalSnapshot snapshot,
                                                           final CanonicalGaussianTransition out) {
        final int[] blockStarts = snapshot.getBlockStarts();
        final int[] blockSizes = snapshot.getBlockSizes();

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

    public void fillTransitionMatrix(final BlockBasisGlobalSnapshot snapshot,
                                     final double dt,
                                     final double[] out) {
        snapshot.ensureCurrent();
        structure.checkMatches(snapshot.getStructure());
        expSolver.computeCompressed(snapshot.getBlockDiagonalCompressed(), dt, expD);
        fillTransitionMatrixFromCompressed(
                expD, out, dimension, snapshot.getBlockStarts(), snapshot.getBlockSizes());
    }

    public void fillTransitionCovariance(final BlockBasisGlobalSnapshot snapshot,
                                         final double dt,
                                         final double[] out) {
        computeTransitionMoments(snapshot, dt);
        System.arraycopy(omegaTilde.data, 0, out, 0, dimension * dimension);
    }

    public void fillTransitionMoments(final BlockBasisGlobalSnapshot snapshot,
                                      final double dt,
                                      final double[] transitionMatrixOut,
                                      final double[] covarianceOut) {
        computeTransitionMoments(snapshot, dt);
        fillTransitionMatrixFromCompressed(
                expD, transitionMatrixOut, dimension, snapshot.getBlockStarts(), snapshot.getBlockSizes());
        System.arraycopy(omegaTilde.data, 0, covarianceOut, 0, dimension * dimension);
    }

    public void fillBranchLengthDerivatives(final BlockBasisGlobalSnapshot snapshot,
                                            final double dt,
                                            final double[] transitionMatrixDerivativeOut,
                                            final double[] covarianceDerivativeOut) {
        snapshot.ensureCurrent();
        structure.checkMatches(snapshot.getStructure());
        final int[] blockStarts = snapshot.getBlockStarts();
        final int[] blockSizes = snapshot.getBlockSizes();
        final double[] blockD = snapshot.getBlockDiagonalCompressed();
        final DenseMatrix64F vTilde = snapshot.getStationaryCovarianceDBasis();

        expSolver.computeCompressed(blockD, dt, expD);
        fillCompressedBranchLengthDerivative(blockD, expD, dExpDdt, dimension, blockStarts, blockSizes);
        fillTransitionMatrixFromCompressed(expD, mFlat, dimension, blockStarts, blockSizes);
        fillTransitionMatrixFromCompressed(dExpDdt, transitionMatrixDerivativeOut,
                dimension, blockStarts, blockSizes);

        BlockDiagonalMatrixOps.multiplyBlockDiagonalLeft(
                dExpDdt, vTilde.data, atildeVTilde.data, dimension, blockStarts, blockSizes);
        BlockDiagonalMatrixOps.multiplyRightBlockDiagonalTranspose(
                atildeVTilde.data, expD, covarianceDerivativeOut, dimension, blockStarts, blockSizes);

        BlockDiagonalMatrixOps.multiplyBlockDiagonalLeft(
                expD, vTilde.data, atildeVTilde.data, dimension, blockStarts, blockSizes);
        BlockDiagonalMatrixOps.multiplyRightBlockDiagonalTranspose(
                atildeVTilde.data, dExpDdt, xxFlat, dimension, blockStarts, blockSizes);
        for (int i = 0; i < dimension * dimension; ++i) {
            covarianceDerivativeOut[i] = -(covarianceDerivativeOut[i] + xxFlat[i]);
        }
        symmetrizeFlat(covarianceDerivativeOut, dimension);
    }

    private void computeTransitionMoments(final BlockBasisGlobalSnapshot snapshot,
                                          final double dt) {
        snapshot.ensureCurrent();
        structure.checkMatches(snapshot.getStructure());
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

    private static void symmetrizeFlat(final double[] data, final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = i + 1; j < dimension; ++j) {
                final double value = 0.5 * (data[i * dimension + j] + data[j * dimension + i]);
                data[i * dimension + j] = value;
                data[j * dimension + i] = value;
            }
        }
    }

    private static void fillTransitionMatrixFromCompressed(final double[] compressed,
                                                           final double[] out,
                                                           final int dimension,
                                                           final int[] blockStarts,
                                                           final int[] blockSizes) {
        if (out.length < dimension * dimension) {
            throw new IllegalArgumentException("Transition matrix buffer is too small.");
        }
        Arrays.fill(out, 0.0);
        int block2x2Index = 0;
        final int upperBase = dimension;
        final int lowerBase = dimension + BlockDiagonalStructure.numTwoByTwoBlocks(blockSizes);
        for (int b = 0; b < blockStarts.length; ++b) {
            final int start = blockStarts[b];
            if (blockSizes[b] == 1) {
                out[start * dimension + start] = compressed[start];
            } else {
                out[start * dimension + start] = compressed[start];
                out[(start + 1) * dimension + start + 1] = compressed[start + 1];
                out[start * dimension + start + 1] = compressed[upperBase + block2x2Index];
                out[(start + 1) * dimension + start] = compressed[lowerBase + block2x2Index];
                block2x2Index++;
            }
        }
    }

    private static void fillCompressedBranchLengthDerivative(final double[] blockD,
                                                             final double[] expD,
                                                             final double[] out,
                                                             final int dimension,
                                                             final int[] blockStarts,
                                                             final int[] blockSizes) {
        Arrays.fill(out, 0.0);
        int block2x2Index = 0;
        final int rawUpperBase = dimension;
        final int rawLowerBase = dimension + dimension - 1;
        final int activeUpperBase = dimension;
        final int activeLowerBase = dimension + BlockDiagonalStructure.numTwoByTwoBlocks(blockSizes);
        for (int b = 0; b < blockStarts.length; ++b) {
            final int start = blockStarts[b];
            if (blockSizes[b] == 1) {
                out[start] = -blockD[start] * expD[start];
            } else {
                final int rawUpperIndex = rawUpperBase + start;
                final int rawLowerIndex = rawLowerBase + start;
                final int activeUpperIndex = activeUpperBase + block2x2Index;
                final int activeLowerIndex = activeLowerBase + block2x2Index;
                final double b00 = blockD[start];
                final double b11 = blockD[start + 1];
                final double b01 = blockD[rawUpperIndex];
                final double b10 = blockD[rawLowerIndex];
                final double e00 = expD[start];
                final double e11 = expD[start + 1];
                final double e01 = expD[activeUpperIndex];
                final double e10 = expD[activeLowerIndex];
                out[start] = -(b00 * e00 + b01 * e10);
                out[start + 1] = -(b10 * e01 + b11 * e11);
                out[activeUpperIndex] = -(b00 * e01 + b01 * e11);
                out[activeLowerIndex] = -(b10 * e00 + b11 * e10);
                block2x2Index++;
            }
        }
    }
}
