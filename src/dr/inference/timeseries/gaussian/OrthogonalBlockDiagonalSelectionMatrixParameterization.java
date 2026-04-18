package dr.inference.timeseries.gaussian;

import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalExpSolver;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalFrechetHelper;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalLyapunovAdjointHelper;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalLyapunovSolver;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.OrthogonalMatrixProvider;
import dr.inference.timeseries.engine.gaussian.CanonicalBranchMessageContribution;
import dr.inference.timeseries.engine.gaussian.CanonicalLocalTransitionAdjoints;
import dr.inference.timeseries.engine.gaussian.CanonicalTransitionAdjointUtils;
import dr.inference.timeseries.representation.CanonicalGaussianTransition;
import dr.inference.timeseries.representation.CanonicalGaussianUtils;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.Arrays;

/**
 * Orthogonal block-diagonal parametrization using block-space helpers for the
 * forward path and the orthogonal native backward path.
 */
public final class OrthogonalBlockDiagonalSelectionMatrixParameterization extends DenseSelectionMatrixParameterization {

    private static final double SPD_JITTER_RELATIVE = 1.0e-14;
    private static final double SPD_JITTER_ABSOLUTE = 1.0e-14;
    private static final String NATIVE_FORCE_DENSE_ADJOINT_EXP_PROPERTY =
            "beast.experimental.nativeForceDenseAdjointExp";
    private static final String TRANSPOSE_NATIVE_FRECHET_INPUT_PROPERTY =
            "beast.experimental.transposeNativeFrechetInput";
    private static final String DEBUG_NATIVE_R_CONSISTENCY_PROPERTY =
            "beast.debug.nativeRConsistency";
    private static final String DEBUG_PD_PIVOT_FLOOR_PROPERTY =
            "beast.debug.pdPivotFloor";

    private final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter;
    private final OrthogonalMatrixProvider orthogonalRotation;
    private final BlockDiagonalExpSolver expSolver;
    private final BlockDiagonalFrechetHelper frechetHelper;
    private final BlockDiagonalLyapunovSolver lyapunovSolver;
    private final BlockDiagonalLyapunovAdjointHelper lyapunovAdjointHelper;
    private final double[] rData;
    private final double[] rtData;
    private final double[] blockDParams;
    private final double[] cachedRData;
    private final double[] cachedRtData;
    private final double[] cachedBlockDParams;
    private final DenseMatrix64F expD;
    private final DenseMatrix64F rMatrix;
    private final DenseMatrix64F rtMatrix;
    private final DenseMatrix64F transitionMatrix;
    private final DenseMatrix64F qMatrix;
    private final DenseMatrix64F qDBasis;
    private final DenseMatrix64F stationaryCovDBasis;
    private final DenseMatrix64F transitionCovDBasis;
    private final DenseMatrix64F transitionCovariance;
    private final double[] transitionMatrixArrayScratch;
    private final double[] transitionCovarianceArrayScratch;
    private final double[] precisionFlat;
    private final double[] transitionOffsetScratch;
    private final double[][] scaledNegativeBlockDScratch;
    private final double[][] denseAdjointScratch;
    private final DenseMatrix64F upstreamF;
    private final DenseMatrix64F upstreamFD;
    private final DenseMatrix64F gradD;
    private final DenseMatrix64F gradR;
    private final DenseMatrix64F gV;
    private final DenseMatrix64F hDBasis;
    private final DenseMatrix64F gS;
    private final DenseMatrix64F yAdjoint;
    private final DenseMatrix64F gECov;
    private final DenseMatrix64F temp1;
    private final DenseMatrix64F temp2;
    private final DenseMatrix64F temp3;
    private final double[][] workMatrix;
    private final double[] choleskyScratch;
    private final double[] lowerInverseScratch;
    private final CanonicalTransitionAdjointUtils.Workspace canonicalAdjointWorkspace;
    private final double[] tempVector1;
    private final double[] tempVector2;
    private boolean expCacheValid;
    private boolean basisCacheValid;
    private double cachedExpDt;
    private double cachedBasisDt;

    public static final class PreparedBranchBasis {
        private final int dimension;
        private double dt;
        private final double[] stationaryMean;
        private final double[] blockDParams;
        private final DenseMatrix64F expD;
        private final DenseMatrix64F rMatrix;
        private final DenseMatrix64F rtMatrix;
        private final DenseMatrix64F transitionMatrix;
        private final double[][] workMatrix;

        private PreparedBranchBasis(final int dimension,
                                    final int blockDParamDimension) {
            this.dimension = dimension;
            this.dt = Double.NaN;
            this.stationaryMean = new double[dimension];
            this.blockDParams = new double[blockDParamDimension];
            this.expD = new DenseMatrix64F(dimension, dimension);
            this.rMatrix = new DenseMatrix64F(dimension, dimension);
            this.rtMatrix = new DenseMatrix64F(dimension, dimension);
            this.transitionMatrix = new DenseMatrix64F(dimension, dimension);
            this.workMatrix = new double[dimension][dimension];
        }
    }

    public static final class BranchGradientWorkspace {
        private final BlockDiagonalFrechetHelper frechetHelper;
        private final BlockDiagonalLyapunovSolver lyapunovSolver;
        private final BlockDiagonalLyapunovAdjointHelper lyapunovAdjointHelper;
        private final DenseMatrix64F qMatrix;
        private final DenseMatrix64F qDBasis;
        private final DenseMatrix64F stationaryCovDBasis;
        private final DenseMatrix64F transitionCovDBasis;
        private final DenseMatrix64F transitionCovariance;
        private final double[] transitionCovarianceArrayScratch;
        private final double[] transitionOffsetScratch;
        private final double[][] scaledNegativeBlockDScratch;
        private final double[][] denseAdjointScratch;
        private final DenseMatrix64F upstreamF;
        private final DenseMatrix64F upstreamFD;
        private final DenseMatrix64F gradD;
        private final DenseMatrix64F gradR;
        private final DenseMatrix64F gV;
        private final DenseMatrix64F hDBasis;
        private final DenseMatrix64F gS;
        private final DenseMatrix64F yAdjoint;
        private final DenseMatrix64F gECov;
        private final DenseMatrix64F temp1;
        private final DenseMatrix64F temp2;
        private final DenseMatrix64F temp3;
        private final double[] choleskyScratch;
        private final double[] lowerInverseScratch;
        private final double[] tempVector1;
        private final double[] tempVector2;

        public BranchGradientWorkspace(final int dimension,
                                       final int[] blockStarts,
                                       final int[] blockSizes) {
            final BlockDiagonalExpSolver.BlockStructure structure =
                    new BlockDiagonalExpSolver.BlockStructure(dimension, blockStarts, blockSizes);
            this.frechetHelper = new BlockDiagonalFrechetHelper(structure);
            this.lyapunovSolver = new BlockDiagonalLyapunovSolver(dimension, blockStarts, blockSizes);
            this.lyapunovAdjointHelper = new BlockDiagonalLyapunovAdjointHelper(dimension, lyapunovSolver);
            this.qMatrix = new DenseMatrix64F(dimension, dimension);
            this.qDBasis = new DenseMatrix64F(dimension, dimension);
            this.stationaryCovDBasis = new DenseMatrix64F(dimension, dimension);
            this.transitionCovDBasis = new DenseMatrix64F(dimension, dimension);
            this.transitionCovariance = new DenseMatrix64F(dimension, dimension);
            this.transitionCovarianceArrayScratch = new double[dimension * dimension];
            this.transitionOffsetScratch = new double[dimension];
            this.scaledNegativeBlockDScratch = new double[dimension][dimension];
            this.denseAdjointScratch = new double[dimension][dimension];
            this.upstreamF = new DenseMatrix64F(dimension, dimension);
            this.upstreamFD = new DenseMatrix64F(dimension, dimension);
            this.gradD = new DenseMatrix64F(dimension, dimension);
            this.gradR = new DenseMatrix64F(dimension, dimension);
            this.gV = new DenseMatrix64F(dimension, dimension);
            this.hDBasis = new DenseMatrix64F(dimension, dimension);
            this.gS = new DenseMatrix64F(dimension, dimension);
            this.yAdjoint = new DenseMatrix64F(dimension, dimension);
            this.gECov = new DenseMatrix64F(dimension, dimension);
            this.temp1 = new DenseMatrix64F(dimension, dimension);
            this.temp2 = new DenseMatrix64F(dimension, dimension);
            this.temp3 = new DenseMatrix64F(dimension, dimension);
            this.choleskyScratch = new double[dimension * dimension];
            this.lowerInverseScratch = new double[dimension * dimension];
            this.tempVector1 = new double[dimension];
            this.tempVector2 = new double[dimension];
        }
    }

    public OrthogonalBlockDiagonalSelectionMatrixParameterization(
            final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
            final OrthogonalMatrixProvider orthogonalRotation) {
        super(blockParameter);
        this.blockParameter = blockParameter;
        this.orthogonalRotation = orthogonalRotation;
        this.expSolver = new BlockDiagonalExpSolver(
                new BlockDiagonalExpSolver.BlockStructure(
                        blockParameter.getRowDimension(),
                        blockParameter.getBlockStarts(),
                        blockParameter.getBlockSizes()));
        this.frechetHelper = new BlockDiagonalFrechetHelper(expSolver.getStructure());
        this.lyapunovSolver = new BlockDiagonalLyapunovSolver(
                getDimension(),
                blockParameter.getBlockStarts(),
                blockParameter.getBlockSizes());
        this.lyapunovAdjointHelper = new BlockDiagonalLyapunovAdjointHelper(getDimension(), lyapunovSolver);

        final int d = getDimension();
        this.rData = new double[d * d];
        this.rtData = new double[d * d];
        this.blockDParams = new double[blockParameter.getTridiagonalDDimension()];
        this.cachedRData = new double[d * d];
        this.cachedRtData = new double[d * d];
        this.cachedBlockDParams = new double[blockParameter.getTridiagonalDDimension()];
        this.expD = new DenseMatrix64F(d, d);
        this.rMatrix = new DenseMatrix64F(d, d);
        this.rtMatrix = new DenseMatrix64F(d, d);
        this.transitionMatrix = new DenseMatrix64F(d, d);
        this.qMatrix = new DenseMatrix64F(d, d);
        this.qDBasis = new DenseMatrix64F(d, d);
        this.stationaryCovDBasis = new DenseMatrix64F(d, d);
        this.transitionCovDBasis = new DenseMatrix64F(d, d);
        this.transitionCovariance = new DenseMatrix64F(d, d);
        this.transitionMatrixArrayScratch = new double[d * d];
        this.transitionCovarianceArrayScratch = new double[d * d];
        this.precisionFlat = new double[d * d];
        this.transitionOffsetScratch = new double[d];
        this.scaledNegativeBlockDScratch = new double[d][d];
        this.denseAdjointScratch = new double[d][d];
        this.upstreamF = new DenseMatrix64F(d, d);
        this.upstreamFD = new DenseMatrix64F(d, d);
        this.gradD = new DenseMatrix64F(d, d);
        this.gradR = new DenseMatrix64F(d, d);
        this.gV = new DenseMatrix64F(d, d);
        this.hDBasis = new DenseMatrix64F(d, d);
        this.gS = new DenseMatrix64F(d, d);
        this.yAdjoint = new DenseMatrix64F(d, d);
        this.gECov = new DenseMatrix64F(d, d);
        this.temp1 = new DenseMatrix64F(d, d);
        this.temp2 = new DenseMatrix64F(d, d);
        this.temp3 = new DenseMatrix64F(d, d);
        this.workMatrix = new double[d][d];
        this.choleskyScratch = new double[d * d];
        this.lowerInverseScratch = new double[d * d];
        this.canonicalAdjointWorkspace = new CanonicalTransitionAdjointUtils.Workspace(d);
        this.tempVector1 = new double[d];
        this.tempVector2 = new double[d];
        this.expCacheValid = false;
        this.basisCacheValid = false;
        this.cachedExpDt = Double.NaN;
        this.cachedBasisDt = Double.NaN;
    }

    @Override
    public void fillTransitionMatrix(final double dt, final double[][] out) {
        refreshBasisCaches(dt);
        copyDenseMatrixToArray(transitionMatrix, out);
    }

    public void fillTransitionCovariance(final MatrixParameterInterface diffusionMatrix,
                                         final double dt,
                                         final double[][] out) {
        fillTransitionCovarianceMatrix(diffusionMatrix, dt, transitionCovariance);
        copyDenseMatrixToArray(transitionCovariance, out);
    }

    public void fillCanonicalTransition(final MatrixParameterInterface diffusionMatrix,
                                        final double[] stationaryMean,
                                        final double dt,
                                        final CanonicalGaussianTransition out) {
        refreshBasisCaches(dt);
        fillTransitionCovarianceMatrix(diffusionMatrix, dt, transitionCovariance);
        fillCanonicalTransitionDirect(stationaryMean, out);
    }

    public PreparedBranchBasis prepareBranchBasis(final double dt,
                                                  final double[] stationaryMean) {
        final PreparedBranchBasis prepared = createPreparedBranchBasis();
        prepareBranchBasis(dt, stationaryMean, prepared);
        return prepared;
    }

    public PreparedBranchBasis createPreparedBranchBasis() {
        final int dimension = getDimension();
        return new PreparedBranchBasis(dimension, blockParameter.getTridiagonalDDimension());
    }

    public void prepareBranchBasis(final double dt,
                                   final double[] stationaryMean,
                                   final PreparedBranchBasis prepared) {
        final int dimension = getDimension();
        if (prepared.dimension != dimension) {
            throw new IllegalArgumentException(
                    "prepared basis dimension must be " + dimension + " but is " + prepared.dimension);
        }
        if (stationaryMean.length != dimension) {
            throw new IllegalArgumentException(
                    "stationaryMean must have length " + dimension + " but has " + stationaryMean.length);
        }

        prepared.dt = dt;
        System.arraycopy(stationaryMean, 0, prepared.stationaryMean, 0, dimension);
        orthogonalRotation.fillOrthogonalMatrix(prepared.rMatrix.data);
        orthogonalRotation.fillOrthogonalTranspose(prepared.rtMatrix.data);
        blockParameter.fillBlockDiagonalElements(prepared.blockDParams);
        expSolver.compute(prepared.blockDParams, dt, prepared.expD);
        multiplyRowMajor(
                prepared.rMatrix.data,
                prepared.expD.data,
                dimension,
                prepared.workMatrix);
        multiply(
                prepared.workMatrix,
                prepared.rtMatrix.data,
                dimension,
                prepared.transitionMatrix);
    }

    public BranchGradientWorkspace createBranchGradientWorkspace() {
        return new BranchGradientWorkspace(
                getDimension(),
                blockParameter.getBlockStarts(),
                blockParameter.getBlockSizes());
    }

    public void fillCanonicalTransitionPrepared(final PreparedBranchBasis prepared,
                                                final MatrixParameterInterface diffusionMatrix,
                                                final BranchGradientWorkspace workspace,
                                                final CanonicalGaussianTransition out) {
        fillTransitionCovarianceMatrixPrepared(prepared, diffusionMatrix, workspace, workspace.transitionCovariance);
        fillCanonicalTransitionDirectPrepared(prepared, workspace.transitionCovariance, workspace, out);
    }

    public void accumulateNativeGradientFromAdjointsPrepared(final PreparedBranchBasis prepared,
                                                             final MatrixParameterInterface diffusionMatrix,
                                                             final CanonicalLocalTransitionAdjoints localAdjoints,
                                                             final BranchGradientWorkspace workspace,
                                                             final double[] compressedDAccumulator,
                                                             final double[][] rotationAccumulator) {
        fillTransitionCovarianceMatrixPrepared(prepared, diffusionMatrix, workspace, workspace.transitionCovariance);

        accumulateNativeGradientFromTransitionPrepared(
                prepared,
                localAdjoints.dLogL_dF,
                localAdjoints.dLogL_df,
                workspace,
                compressedDAccumulator,
                rotationAccumulator);
        if (!isFinite(compressedDAccumulator) || !isFinite(rotationAccumulator)) {
            throw new IllegalStateException(
                    "Non-finite orthogonal native transition contribution at dt=" + prepared.dt);
        }

        accumulateNativeGradientFromCovarianceStationaryPrepared(
                prepared,
                localAdjoints.dLogL_dOmega,
                workspace,
                compressedDAccumulator,
                rotationAccumulator);
        if (!isFinite(compressedDAccumulator) || !isFinite(rotationAccumulator)) {
            throw new IllegalStateException(
                    "Non-finite orthogonal native covariance contribution at dt=" + prepared.dt);
        }
    }

    public void accumulateMeanGradientPrepared(final PreparedBranchBasis prepared,
                                               final double[] dLogL_df,
                                               final double[] gradientAccumulator,
                                               final BranchGradientWorkspace workspace) {
        multiplyRowMajorVector(prepared.rtMatrix.data, dLogL_df, workspace.tempVector1, prepared.dimension);
        applyBlockDiagonalTransposeExpPrepared(prepared, workspace.tempVector1, workspace.tempVector2);
        multiplyRowMajorVector(prepared.rMatrix.data, workspace.tempVector2, workspace.tempVector1, prepared.dimension);
        if (gradientAccumulator.length == 1) {
            double sum = 0.0;
            for (int i = 0; i < prepared.dimension; ++i) {
                sum += dLogL_df[i] - workspace.tempVector1[i];
            }
            gradientAccumulator[0] += sum;
            return;
        }
        if (gradientAccumulator.length != prepared.dimension) {
            throw new IllegalArgumentException(
                    "Stationary-mean gradient length must be 1 or "
                            + prepared.dimension + ", found " + gradientAccumulator.length);
        }
        for (int i = 0; i < prepared.dimension; ++i) {
            gradientAccumulator[i] += dLogL_df[i] - workspace.tempVector1[i];
        }
    }

    public void accumulateDiffusionGradientPrepared(final PreparedBranchBasis prepared,
                                                    final double[][] dLogL_dV,
                                                    final double[] gradientAccumulator,
                                                    final BranchGradientWorkspace workspace) {
        fillSymmetricDenseMatrix(dLogL_dV, workspace.gV);

        CommonOps.mult(prepared.rtMatrix, workspace.gV, workspace.temp1);
        CommonOps.mult(workspace.temp1, prepared.rMatrix, workspace.hDBasis);
        symmetrize(workspace.hDBasis);

        CommonOps.multTransA(prepared.expD, workspace.hDBasis, workspace.temp1);
        CommonOps.mult(workspace.temp1, prepared.expD, workspace.gS);
        CommonOps.subtract(workspace.hDBasis, workspace.gS, workspace.gS);
        symmetrize(workspace.gS);

        workspace.yAdjoint.set(workspace.gS);
        CommonOps.scale(-1.0, workspace.yAdjoint);
        workspace.lyapunovAdjointHelper.solveAdjointInDBasis(
                workspace.yAdjoint, prepared.blockDParams, workspace.yAdjoint);

        CommonOps.mult(prepared.rMatrix, workspace.yAdjoint, workspace.temp1);
        CommonOps.mult(workspace.temp1, prepared.rtMatrix, workspace.temp2);
        symmetrize(workspace.temp2);

        final int dimension = prepared.dimension;
        final double[] data = workspace.temp2.data;
        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                gradientAccumulator[rowOffset + j] += data[rowOffset + j];
            }
        }
    }

    public void fillCanonicalLocalAdjoints(final MatrixParameterInterface diffusionMatrix,
                                           final double[] stationaryMean,
                                           final double dt,
                                           final CanonicalBranchMessageContribution contribution,
                                           final CanonicalLocalTransitionAdjoints out) {
        refreshBasisCaches(dt);
        fillTransitionCovarianceMatrix(diffusionMatrix, dt, transitionCovariance);
        copyDenseMatrixToFlat(transitionMatrix, transitionMatrixArrayScratch);
        fillTransitionOffset(stationaryMean, transitionOffsetScratch);

        copyAndInvertPositiveDefiniteFlat(transitionCovariance, transitionCovarianceArrayScratch, precisionFlat);
        CanonicalTransitionAdjointUtils.fillFromMoments(
                precisionFlat,
                transitionCovarianceArrayScratch,
                transitionMatrixArrayScratch,
                transitionOffsetScratch,
                contribution,
                canonicalAdjointWorkspace,
                out);
    }

    private void fillTransitionCovarianceMatrix(final MatrixParameterInterface diffusionMatrix,
                                                final double dt,
                                                final DenseMatrix64F out) {
        refreshBasisCaches(dt);
        fillDenseMatrix(diffusionMatrix, qMatrix);
        CommonOps.mult(rtMatrix, qMatrix, temp1);
        CommonOps.mult(temp1, rMatrix, qDBasis);
        lyapunovSolver.solve(blockDParams, qDBasis, stationaryCovDBasis);
        CommonOps.mult(expD, stationaryCovDBasis, temp1);
        CommonOps.multTransB(temp1, expD, transitionCovDBasis);
        CommonOps.subtract(stationaryCovDBasis, transitionCovDBasis, transitionCovDBasis);
        CommonOps.mult(rMatrix, transitionCovDBasis, temp1);
        CommonOps.mult(temp1, rtMatrix, out);
        symmetrize(out);
    }

    private void fillCanonicalTransitionDirect(final double[] stationaryMean,
                                               final CanonicalGaussianTransition out) {
        final int d = getDimension();
        final double[] transitionData = transitionMatrix.data;
        final double[] transitionOffset = transitionOffsetScratch;
        fillTransitionOffset(stationaryMean, transitionOffset);

        final double logDet = copyAndInvertPositiveDefiniteFlat(
                transitionCovariance,
                transitionCovarianceArrayScratch,
                out.precisionYY);

        // J_yx = -P F
        for (int i = 0; i < d; ++i) {
            final int iOff = i * d;
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += out.precisionYY[iOff + k] * transitionData[k * d + j];
                }
                out.precisionYX[iOff + j] = -sum;
            }
        }

        // J_xy = J_yx^T
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                out.precisionXY[i * d + j] = out.precisionYX[j * d + i];
            }
        }

        // J_xx = F^T P F = -(F^T J_yx)
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum -= transitionData[k * d + i] * out.precisionYX[k * d + j];
                }
                out.precisionXX[i * d + j] = sum;
            }
        }

        // h_y = P f
        for (int i = 0; i < d; ++i) {
            double sum = 0.0;
            final int iOff = i * d;
            for (int j = 0; j < d; ++j) {
                sum += out.precisionYY[iOff + j] * transitionOffset[j];
            }
            out.informationY[i] = sum;
        }

        // h_x = -F^T h_y
        for (int i = 0; i < d; ++i) {
            double sum = 0.0;
            for (int j = 0; j < d; ++j) {
                sum += transitionData[j * d + i] * out.informationY[j];
            }
            out.informationX[i] = -sum;
        }

        double quadratic = 0.0;
        for (int i = 0; i < d; ++i) {
            quadratic += transitionOffset[i] * out.informationY[i];
        }
        out.logNormalizer = 0.5 * (d * Math.log(2.0 * Math.PI) + logDet + quadratic);
    }

    public void accumulateNativeGradientFromTransition(final double dt,
                                                       final double[] stationaryMean,
                                                       final double[][] dLogL_dF,
                                                       final double[] dLogL_df,
                                                       final double[] compressedDAccumulator,
                                                       final double[][] rotationAccumulator) {
        copySquareMatrixToFlat(dLogL_dF, transitionMatrixArrayScratch);
        accumulateNativeGradientFromTransition(
                dt, stationaryMean, transitionMatrixArrayScratch, dLogL_df,
                compressedDAccumulator, rotationAccumulator);
    }

    public void accumulateNativeGradientFromTransition(final double dt,
                                                       final double[] stationaryMean,
                                                       final double[] dLogL_dF,
                                                       final double[] dLogL_df,
                                                       final double[] compressedDAccumulator,
                                                       final double[][] rotationAccumulator) {
        refreshBasisCaches(dt);
        accumulateNativeGradientFromTransitionCached(
                dt, stationaryMean, dLogL_dF, dLogL_df, compressedDAccumulator, rotationAccumulator);
    }

    public void accumulateNativeGradientFromCovarianceStationary(final MatrixParameterInterface diffusionMatrix,
                                                                 final double dt,
                                                                 final double[][] dLogL_dV,
                                                                 final double[] compressedDAccumulator,
                                                                 final double[][] rotationAccumulator) {
        copySquareMatrixToFlat(dLogL_dV, transitionCovarianceArrayScratch);
        accumulateNativeGradientFromCovarianceStationary(
                diffusionMatrix, dt, transitionCovarianceArrayScratch,
                compressedDAccumulator, rotationAccumulator);
    }

    public void accumulateNativeGradientFromCovarianceStationary(final MatrixParameterInterface diffusionMatrix,
                                                                 final double dt,
                                                                 final double[] dLogL_dV,
                                                                 final double[] compressedDAccumulator,
                                                                 final double[][] rotationAccumulator) {
        refreshBasisCaches(dt);
        fillDenseMatrix(diffusionMatrix, qMatrix);
        CommonOps.mult(rtMatrix, qMatrix, temp1);
        CommonOps.mult(temp1, rMatrix, qDBasis);
        lyapunovSolver.solve(blockDParams, qDBasis, stationaryCovDBasis);
        CommonOps.mult(expD, stationaryCovDBasis, temp1);
        CommonOps.multTransB(temp1, expD, transitionCovDBasis);
        CommonOps.subtract(stationaryCovDBasis, transitionCovDBasis, transitionCovDBasis);
        symmetrize(transitionCovDBasis);
        accumulateNativeGradientFromCovarianceStationaryCached(
                dt, dLogL_dV, compressedDAccumulator, rotationAccumulator);
    }

    public void accumulateDiffusionGradient(final MatrixParameterInterface diffusionMatrix,
                                            final double dt,
                                            final double[][] dLogL_dV,
                                            final double[] gradientAccumulator) {
        refreshBasisCaches(dt);
        fillSymmetricDenseMatrix(dLogL_dV, gV);

        // Map covariance adjoint to block-D basis: G_V^D = R^T G_V R.
        CommonOps.mult(rtMatrix, gV, temp1);
        CommonOps.mult(temp1, rMatrix, hDBasis);
        symmetrize(hDBasis);

        // For V = R (S - E S E^T) R^T with D S + S D^T = Q_D, the exact adjoint is:
        //   G_S = G_V^D - E^T G_V^D E
        //   D^T Y + Y D = G_S
        //   G_Q^D = Y
        //   G_Q = R G_Q^D R^T
        CommonOps.multTransA(expD, hDBasis, temp1); // E^T G_V^D
        CommonOps.mult(temp1, expD, gS);            // E^T G_V^D E
        CommonOps.subtract(hDBasis, gS, gS);        // G_S
        symmetrize(gS);

        // Helper solves D^T Y + Y D = -hBlock, so pass hBlock = -G_S.
        yAdjoint.set(gS);
        CommonOps.scale(-1.0, yAdjoint);
        lyapunovAdjointHelper.solveAdjointInDBasis(yAdjoint, blockDParams, yAdjoint);

        CommonOps.mult(rMatrix, yAdjoint, temp1);
        CommonOps.mult(temp1, rtMatrix, temp2);
        symmetrize(temp2);

        final int d = getDimension();
        final double[] data = temp2.data;
        for (int i = 0; i < d; ++i) {
            final int rowOffset = i * d;
            for (int j = 0; j < d; ++j) {
                gradientAccumulator[rowOffset + j] += data[rowOffset + j];
            }
        }
    }

    public void accumulateNativeGradientFromCanonicalContribution(final MatrixParameterInterface diffusionMatrix,
                                                                  final double[] stationaryMean,
                                                                  final double dt,
                                                                  final CanonicalBranchMessageContribution contribution,
                                                                  final CanonicalLocalTransitionAdjoints localAdjoints,
                                                                  final double[] compressedDAccumulator,
                                                                  final double[][] rotationAccumulator) {
        refreshBasisCaches(dt);
        fillTransitionCovarianceMatrix(diffusionMatrix, dt, transitionCovariance);
        copyDenseMatrixToFlat(transitionMatrix, transitionMatrixArrayScratch);
        fillTransitionOffset(stationaryMean, transitionOffsetScratch);

        copyAndInvertPositiveDefiniteFlat(transitionCovariance, transitionCovarianceArrayScratch, precisionFlat);
        CanonicalTransitionAdjointUtils.fillFromMoments(
                precisionFlat,
                transitionCovarianceArrayScratch,
                transitionMatrixArrayScratch,
                transitionOffsetScratch,
                contribution,
                canonicalAdjointWorkspace,
                localAdjoints);

        accumulateNativeGradientFromTransitionCached(
                dt,
                stationaryMean,
                localAdjoints.dLogL_dF,
                localAdjoints.dLogL_df,
                compressedDAccumulator,
                rotationAccumulator);
        if (!isFinite(compressedDAccumulator) || !isFinite(rotationAccumulator)) {
            throw new IllegalStateException("Non-finite orthogonal native transition contribution at dt=" + dt);
        }
        accumulateNativeGradientFromCovarianceStationaryCached(
                dt,
                localAdjoints.dLogL_dOmega,
                compressedDAccumulator,
                rotationAccumulator);
        if (!isFinite(compressedDAccumulator) || !isFinite(rotationAccumulator)) {
            throw new IllegalStateException("Non-finite orthogonal native covariance contribution at dt=" + dt);
        }
    }

    public void accumulateNativeGradientFromAdjoints(final MatrixParameterInterface diffusionMatrix,
                                                     final double[] stationaryMean,
                                                     final double dt,
                                                     final CanonicalLocalTransitionAdjoints localAdjoints,
                                                     final double[] compressedDAccumulator,
                                                     final double[][] rotationAccumulator) {
        refreshBasisCaches(dt);
        fillTransitionCovarianceMatrix(diffusionMatrix, dt, transitionCovariance);

        accumulateNativeGradientFromTransitionCached(
                dt,
                stationaryMean,
                localAdjoints.dLogL_dF,
                localAdjoints.dLogL_df,
                compressedDAccumulator,
                rotationAccumulator);
        if (!isFinite(compressedDAccumulator) || !isFinite(rotationAccumulator)) {
            throw new IllegalStateException("Non-finite orthogonal native transition contribution at dt=" + dt);
        }

        accumulateNativeGradientFromCovarianceStationaryCached(
                dt,
                localAdjoints.dLogL_dOmega,
                compressedDAccumulator,
                rotationAccumulator);
        if (!isFinite(compressedDAccumulator) || !isFinite(rotationAccumulator)) {
            throw new IllegalStateException("Non-finite orthogonal native covariance contribution at dt=" + dt);
        }
    }

    public void accumulateMeanGradient(final double dt,
                                       final double[] dLogL_df,
                                       final double[] gradientAccumulator) {
        refreshMeanGradientCaches(dt);
        orthogonalRotation.applyOrthogonalTranspose(dLogL_df, tempVector1);
        applyBlockDiagonalTransposeExp(tempVector1, tempVector2);
        orthogonalRotation.applyOrthogonal(tempVector2, tempVector1);
        for (int i = 0; i < getDimension(); ++i) {
            gradientAccumulator[i] += dLogL_df[i] - tempVector1[i];
        }
    }

    public double accumulateScalarMeanGradient(final double dt,
                                               final double[] dLogL_df) {
        refreshMeanGradientCaches(dt);
        orthogonalRotation.applyOrthogonalTranspose(dLogL_df, tempVector1);
        applyBlockDiagonalTransposeExp(tempVector1, tempVector2);
        orthogonalRotation.applyOrthogonal(tempVector2, tempVector1);
        double sum = 0.0;
        for (int i = 0; i < getDimension(); ++i) {
            sum += dLogL_df[i] - tempVector1[i];
        }
        return sum;
    }

    private void accumulateNativeGradientFromTransitionCached(final double dt,
                                                              final double[] stationaryMean,
                                                              final double[] dLogL_dF,
                                                              final double[] dLogL_df,
                                                              final double[] compressedDAccumulator,
                                                              final double[][] rotationAccumulator) {
        fillTotalUpstreamOnTransition(stationaryMean, dLogL_dF, dLogL_df, upstreamF);
        final boolean forceDenseAdjointExp = Boolean.getBoolean(NATIVE_FORCE_DENSE_ADJOINT_EXP_PROPERTY);

        CommonOps.mult(rtMatrix, upstreamF, temp1);
        CommonOps.mult(temp1, rMatrix, upstreamFD);
        if (forceDenseAdjointExp) {
            fillScaledNegativeBlockDMatrix(dt, scaledNegativeBlockDScratch);
            copyDenseMatrixToArray(upstreamFD, denseAdjointScratch);
            MatrixExponentialUtils.adjointExp(
                    scaledNegativeBlockDScratch,
                    denseAdjointScratch,
                    denseAdjointScratch);
            fillDenseMatrix(denseAdjointScratch, gradD);
            transposeInPlace(gradD);
            CommonOps.scale(-dt, gradD);
        } else {
            final DenseMatrix64F frechetInputTransition;
            if (Boolean.getBoolean(TRANSPOSE_NATIVE_FRECHET_INPUT_PROPERTY)) {
                CommonOps.transpose(upstreamFD, temp3);
                frechetInputTransition = temp3;
            } else {
                frechetInputTransition = upstreamFD;
            }
            frechetHelper.frechetAdjointExpInDBasis(blockDParams, frechetInputTransition, dt, gradD);
            transposeInPlace(gradD);
            if (!isFinite(gradD.data)) {
                fillScaledNegativeBlockDMatrix(dt, scaledNegativeBlockDScratch);
                copyDenseMatrixToArray(upstreamFD, denseAdjointScratch);
                MatrixExponentialUtils.adjointExp(
                        scaledNegativeBlockDScratch,
                        denseAdjointScratch,
                        denseAdjointScratch);
                fillDenseMatrix(denseAdjointScratch, gradD);
                transposeInPlace(gradD);
                CommonOps.scale(-dt, gradD);
            }
        }
        accumulateCompressedGradient(gradD, compressedDAccumulator);

        CommonOps.multTransB(rMatrix, expD, temp1);     // R * E^T
        CommonOps.mult(upstreamF, temp1, gradR);        // U * R * E^T
        CommonOps.multTransA(upstreamF, rMatrix, temp1);// U^T * R
        CommonOps.mult(temp1, expD, temp2);             // U^T * R * E
        CommonOps.addEquals(gradR, temp2);
        addDenseMatrixToArray(gradR, rotationAccumulator);
    }

    private void fillTransitionCovarianceMatrixPrepared(final PreparedBranchBasis prepared,
                                                        final MatrixParameterInterface diffusionMatrix,
                                                        final BranchGradientWorkspace workspace,
                                                        final DenseMatrix64F out) {
        fillDenseMatrix(diffusionMatrix, workspace.qMatrix);
        CommonOps.mult(prepared.rtMatrix, workspace.qMatrix, workspace.temp1);
        CommonOps.mult(workspace.temp1, prepared.rMatrix, workspace.qDBasis);
        workspace.lyapunovSolver.solve(prepared.blockDParams, workspace.qDBasis, workspace.stationaryCovDBasis);
        CommonOps.mult(prepared.expD, workspace.stationaryCovDBasis, workspace.temp1);
        CommonOps.multTransB(workspace.temp1, prepared.expD, workspace.transitionCovDBasis);
        CommonOps.subtract(workspace.stationaryCovDBasis, workspace.transitionCovDBasis, workspace.transitionCovDBasis);
        symmetrize(workspace.transitionCovDBasis);
        CommonOps.mult(prepared.rMatrix, workspace.transitionCovDBasis, workspace.temp1);
        CommonOps.mult(workspace.temp1, prepared.rtMatrix, out);
        symmetrize(out);
    }

    private void fillCanonicalTransitionDirectPrepared(final PreparedBranchBasis prepared,
                                                       final DenseMatrix64F transitionCovariance,
                                                       final BranchGradientWorkspace workspace,
                                                       final CanonicalGaussianTransition out) {
        final int dimension = prepared.dimension;
        final double[] transitionData = prepared.transitionMatrix.data;
        fillTransitionOffset(prepared.transitionMatrix, prepared.stationaryMean, workspace.transitionOffsetScratch);

        final double logDet = copyAndInvertPositiveDefiniteFlat(
                transitionCovariance,
                workspace.transitionCovarianceArrayScratch,
                out.precisionYY,
                workspace.choleskyScratch,
                workspace.lowerInverseScratch);

        for (int i = 0; i < dimension; ++i) {
            final int iOff = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += out.precisionYY[iOff + k] * transitionData[k * dimension + j];
                }
                out.precisionYX[iOff + j] = -sum;
            }
        }

        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                out.precisionXY[i * dimension + j] = out.precisionYX[j * dimension + i];
            }
        }

        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum -= transitionData[k * dimension + i] * out.precisionYX[k * dimension + j];
                }
                out.precisionXX[i * dimension + j] = sum;
            }
        }

        for (int i = 0; i < dimension; ++i) {
            double sum = 0.0;
            final int iOff = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                sum += out.precisionYY[iOff + j] * workspace.transitionOffsetScratch[j];
            }
            out.informationY[i] = sum;
        }

        for (int i = 0; i < dimension; ++i) {
            double sum = 0.0;
            for (int j = 0; j < dimension; ++j) {
                sum += transitionData[j * dimension + i] * out.informationY[j];
            }
            out.informationX[i] = -sum;
        }

        double quadratic = 0.0;
        for (int i = 0; i < dimension; ++i) {
            quadratic += workspace.transitionOffsetScratch[i] * out.informationY[i];
        }
        out.logNormalizer = 0.5 * (dimension * Math.log(2.0 * Math.PI) + logDet + quadratic);
    }

    private void accumulateNativeGradientFromTransitionPrepared(final PreparedBranchBasis prepared,
                                                                final double[] dLogL_dF,
                                                                final double[] dLogL_df,
                                                                final BranchGradientWorkspace workspace,
                                                                final double[] compressedDAccumulator,
                                                                final double[][] rotationAccumulator) {
        fillTotalUpstreamOnTransition(prepared.stationaryMean, dLogL_dF, dLogL_df, workspace.upstreamF);
        final boolean forceDenseAdjointExp = Boolean.getBoolean(NATIVE_FORCE_DENSE_ADJOINT_EXP_PROPERTY);

        CommonOps.mult(prepared.rtMatrix, workspace.upstreamF, workspace.temp1);
        CommonOps.mult(workspace.temp1, prepared.rMatrix, workspace.upstreamFD);
        if (forceDenseAdjointExp) {
            fillScaledNegativeBlockDMatrix(
                    prepared.blockDParams,
                    prepared.dt,
                    workspace.scaledNegativeBlockDScratch,
                    prepared.dimension);
            copyDenseMatrixToArray(workspace.upstreamFD, workspace.denseAdjointScratch);
            MatrixExponentialUtils.adjointExp(
                    workspace.scaledNegativeBlockDScratch,
                    workspace.denseAdjointScratch,
                    workspace.denseAdjointScratch);
            fillDenseMatrix(workspace.denseAdjointScratch, workspace.gradD);
            transposeInPlace(workspace.gradD);
            CommonOps.scale(-prepared.dt, workspace.gradD);
        } else {
            final DenseMatrix64F frechetInputTransition;
            if (Boolean.getBoolean(TRANSPOSE_NATIVE_FRECHET_INPUT_PROPERTY)) {
                CommonOps.transpose(workspace.upstreamFD, workspace.temp3);
                frechetInputTransition = workspace.temp3;
            } else {
                frechetInputTransition = workspace.upstreamFD;
            }
            workspace.frechetHelper.frechetAdjointExpInDBasis(
                    prepared.blockDParams, frechetInputTransition, prepared.dt, workspace.gradD);
            transposeInPlace(workspace.gradD);
            if (!isFinite(workspace.gradD.data)) {
                fillScaledNegativeBlockDMatrix(
                        prepared.blockDParams,
                        prepared.dt,
                        workspace.scaledNegativeBlockDScratch,
                        prepared.dimension);
                copyDenseMatrixToArray(workspace.upstreamFD, workspace.denseAdjointScratch);
                MatrixExponentialUtils.adjointExp(
                        workspace.scaledNegativeBlockDScratch,
                        workspace.denseAdjointScratch,
                        workspace.denseAdjointScratch);
                fillDenseMatrix(workspace.denseAdjointScratch, workspace.gradD);
                transposeInPlace(workspace.gradD);
                CommonOps.scale(-prepared.dt, workspace.gradD);
            }
        }
        accumulateCompressedGradient(workspace.gradD, compressedDAccumulator);

        CommonOps.multTransB(prepared.rMatrix, prepared.expD, workspace.temp1);
        CommonOps.mult(workspace.upstreamF, workspace.temp1, workspace.gradR);
        CommonOps.multTransA(workspace.upstreamF, prepared.rMatrix, workspace.temp1);
        CommonOps.mult(workspace.temp1, prepared.expD, workspace.temp2);
        CommonOps.addEquals(workspace.gradR, workspace.temp2);
        addDenseMatrixToArray(workspace.gradR, rotationAccumulator);
    }

    private void accumulateNativeGradientFromCovarianceStationaryPrepared(final PreparedBranchBasis prepared,
                                                                          final double[] dLogL_dV,
                                                                          final BranchGradientWorkspace workspace,
                                                                          final double[] compressedDAccumulator,
                                                                          final double[][] rotationAccumulator) {
        final boolean forceDenseAdjointExp = Boolean.getBoolean(NATIVE_FORCE_DENSE_ADJOINT_EXP_PROPERTY);
        fillSymmetricDenseMatrixFlat(dLogL_dV, workspace.gV);

        CommonOps.mult(prepared.rtMatrix, workspace.gV, workspace.temp1);
        CommonOps.mult(workspace.temp1, prepared.rMatrix, workspace.hDBasis);
        symmetrize(workspace.hDBasis);

        CommonOps.multTransA(prepared.expD, workspace.hDBasis, workspace.temp1);
        CommonOps.mult(workspace.temp1, prepared.expD, workspace.gS);
        CommonOps.changeSign(workspace.gS);
        CommonOps.addEquals(workspace.gS, workspace.hDBasis);

        CommonOps.fill(workspace.gradD, 0.0);
        workspace.lyapunovAdjointHelper.accumulateLyapunovContributionInDBasis(
                workspace.stationaryCovDBasis, workspace.gS, prepared.blockDParams, workspace.gradD);
        accumulateCompressedGradient(workspace.gradD, compressedDAccumulator);

        workspace.lyapunovAdjointHelper.solveAdjointInDBasis(
                workspace.gS, prepared.blockDParams, workspace.yAdjoint);

        CommonOps.mult(workspace.qMatrix, prepared.rMatrix, workspace.temp1);
        CommonOps.multTransB(workspace.temp1, workspace.yAdjoint, workspace.gradR);
        CommonOps.multTransA(workspace.qMatrix, prepared.rMatrix, workspace.temp1);
        CommonOps.mult(workspace.temp1, workspace.yAdjoint, workspace.temp2);
        CommonOps.addEquals(workspace.gradR, workspace.temp2);
        CommonOps.scale(-1.0, workspace.gradR);
        addDenseMatrixToArray(workspace.gradR, rotationAccumulator);

        CommonOps.mult(workspace.hDBasis, prepared.expD, workspace.temp1);
        CommonOps.mult(workspace.temp1, workspace.stationaryCovDBasis, workspace.gECov);
        CommonOps.scale(-2.0, workspace.gECov);
        if (forceDenseAdjointExp) {
            fillScaledNegativeBlockDMatrix(
                    prepared.blockDParams,
                    prepared.dt,
                    workspace.scaledNegativeBlockDScratch,
                    prepared.dimension);
            copyDenseMatrixToArray(workspace.gECov, workspace.denseAdjointScratch);
            MatrixExponentialUtils.adjointExp(
                    workspace.scaledNegativeBlockDScratch,
                    workspace.denseAdjointScratch,
                    workspace.denseAdjointScratch);
            fillDenseMatrix(workspace.denseAdjointScratch, workspace.gradD);
            transposeInPlace(workspace.gradD);
            CommonOps.scale(-prepared.dt, workspace.gradD);
        } else {
            final DenseMatrix64F frechetInputCov;
            if (Boolean.getBoolean(TRANSPOSE_NATIVE_FRECHET_INPUT_PROPERTY)) {
                CommonOps.transpose(workspace.gECov, workspace.temp3);
                frechetInputCov = workspace.temp3;
            } else {
                frechetInputCov = workspace.gECov;
            }
            workspace.frechetHelper.frechetAdjointExpInDBasis(
                    prepared.blockDParams, frechetInputCov, prepared.dt, workspace.gradD);
            transposeInPlace(workspace.gradD);
        }
        accumulateCompressedGradient(workspace.gradD, compressedDAccumulator);

        CommonOps.mult(workspace.gV, prepared.rMatrix, workspace.temp1);
        CommonOps.mult(workspace.temp1, workspace.transitionCovDBasis, workspace.gradR);
        CommonOps.multTransA(workspace.gV, prepared.rMatrix, workspace.temp2);
        CommonOps.mult(workspace.temp2, workspace.transitionCovDBasis, workspace.temp3);
        CommonOps.addEquals(workspace.gradR, workspace.temp3);
        addDenseMatrixToArray(workspace.gradR, rotationAccumulator);
    }

    private void accumulateNativeGradientFromCovarianceStationaryCached(final double dt,
                                                                        final double[] dLogL_dV,
                                                                        final double[] compressedDAccumulator,
                                                                        final double[][] rotationAccumulator) {
        final boolean forceDenseAdjointExp = Boolean.getBoolean(NATIVE_FORCE_DENSE_ADJOINT_EXP_PROPERTY);
        fillSymmetricDenseMatrixFlat(dLogL_dV, gV);

        CommonOps.mult(rtMatrix, gV, temp1);
        CommonOps.mult(temp1, rMatrix, hDBasis);
        symmetrize(hDBasis);

        CommonOps.multTransA(expD, hDBasis, temp1);     // E^T H
        CommonOps.mult(temp1, expD, gS);                // E^T H E
        CommonOps.changeSign(gS);
        CommonOps.addEquals(gS, hDBasis);               // H - E^T H E

        CommonOps.fill(gradD, 0.0);
        lyapunovAdjointHelper.accumulateLyapunovContributionInDBasis(
                stationaryCovDBasis, gS, blockDParams, gradD);
        accumulateCompressedGradient(gradD, compressedDAccumulator);

        lyapunovAdjointHelper.solveAdjointInDBasis(gS, blockDParams, yAdjoint);

        // Q-basis contribution: helper returns yAdjoint = -Y, where Y solves
        // D^T Y + Y D = G_S.  Since Q_D = R^T Q R, the raw R gradient is
        // Q R Y^T + Q^T R Y = -(Q R yAdjoint^T + Q^T R yAdjoint).
        CommonOps.mult(qMatrix, rMatrix, temp1);
        CommonOps.multTransB(temp1, yAdjoint, gradR);   // Q R y^T
        CommonOps.multTransA(qMatrix, rMatrix, temp1);  // Q^T R
        CommonOps.mult(temp1, yAdjoint, temp2);         // Q^T R y
        CommonOps.addEquals(gradR, temp2);
        CommonOps.scale(-1.0, gradR);
        addDenseMatrixToArray(gradR, rotationAccumulator);

        // E-inside-V contribution.
        CommonOps.mult(hDBasis, expD, temp1);
        CommonOps.mult(temp1, stationaryCovDBasis, gECov);
        CommonOps.scale(-2.0, gECov);
        if (forceDenseAdjointExp) {
            fillScaledNegativeBlockDMatrix(dt, scaledNegativeBlockDScratch);
            copyDenseMatrixToArray(gECov, denseAdjointScratch);
            MatrixExponentialUtils.adjointExp(
                    scaledNegativeBlockDScratch,
                    denseAdjointScratch,
                    denseAdjointScratch);
            fillDenseMatrix(denseAdjointScratch, gradD);
            transposeInPlace(gradD);
            CommonOps.scale(-dt, gradD);
        } else {
            final DenseMatrix64F frechetInputCov;
            if (Boolean.getBoolean(TRANSPOSE_NATIVE_FRECHET_INPUT_PROPERTY)) {
                CommonOps.transpose(gECov, temp3);
                frechetInputCov = temp3;
            } else {
                frechetInputCov = gECov;
            }
            frechetHelper.frechetAdjointExpInDBasis(blockDParams, frechetInputCov, dt, gradD);
            transposeInPlace(gradD);
        }
        accumulateCompressedGradient(gradD, compressedDAccumulator);

        // Outer basis reconstruction V = R V_D R^T.
        CommonOps.mult(gV, rMatrix, temp1);
        CommonOps.mult(temp1, transitionCovDBasis, gradR);
        CommonOps.multTransA(gV, rMatrix, temp2);
        CommonOps.mult(temp2, transitionCovDBasis, temp3);
        CommonOps.addEquals(gradR, temp3);
        addDenseMatrixToArray(gradR, rotationAccumulator);
    }

    private void fillTransitionOffset(final double[] stationaryMean,
                                      final double[] out) {
        final int dimension = getDimension();
        final double[] transitionData = transitionMatrix.data;
        for (int i = 0; i < dimension; ++i) {
            double transformedMean = 0.0;
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                transformedMean += transitionData[rowOffset + j] * stationaryMean[j];
            }
            out[i] = stationaryMean[i] - transformedMean;
        }
    }

    private static void fillTransitionOffset(final DenseMatrix64F transitionMatrix,
                                             final double[] stationaryMean,
                                             final double[] out) {
        final int dimension = transitionMatrix.numRows;
        final double[] transitionData = transitionMatrix.data;
        for (int i = 0; i < dimension; ++i) {
            double transformedMean = 0.0;
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                transformedMean += transitionData[rowOffset + j] * stationaryMean[j];
            }
            out[i] = stationaryMean[i] - transformedMean;
        }
    }

    private void refreshBasisCaches(final double dt) {
        final int d = getDimension();
        orthogonalRotation.fillOrthogonalMatrix(rData);
        orthogonalRotation.fillOrthogonalTranspose(rtData);
        blockParameter.fillBlockDiagonalElements(blockDParams);

        final boolean expNeedsRefresh = !expCacheValid
                || Double.doubleToLongBits(dt) != Double.doubleToLongBits(cachedExpDt)
                || !Arrays.equals(blockDParams, cachedBlockDParams);
        if (Boolean.getBoolean(DEBUG_NATIVE_R_CONSISTENCY_PROPERTY)) {
            final double[] blockR = new double[d * d];
            final double[] blockRinv = new double[d * d];
            blockParameter.fillRAndRinv(blockR, blockRinv);
            double maxR = 0.0;
            double maxRt = 0.0;
            double maxRinvVsTranspose = 0.0;
            for (int i = 0; i < d * d; ++i) {
                maxR = Math.max(maxR, Math.abs(rData[i] - blockR[i]));
                maxRt = Math.max(maxRt, Math.abs(rtData[i] - blockRinv[i]));
            }
            for (int r = 0; r < d; ++r) {
                for (int c = 0; c < d; ++c) {
                    final int rc = r * d + c;
                    final int cr = c * d + r;
                    maxRinvVsTranspose = Math.max(
                            maxRinvVsTranspose,
                            Math.abs(blockRinv[rc] - blockR[cr]));
                }
            }
            System.err.println("nativeRConsistencyDebug dt=" + dt
                    + " maxAbsRDiff=" + maxR
                    + " maxAbsRtDiff=" + maxRt
                    + " maxAbsRinvVsTranspose=" + maxRinvVsTranspose);
        }
        final boolean basisNeedsRefresh = expNeedsRefresh
                || !basisCacheValid
                || Double.doubleToLongBits(dt) != Double.doubleToLongBits(cachedBasisDt)
                || !Arrays.equals(rData, cachedRData)
                || !Arrays.equals(rtData, cachedRtData);
        if (!basisNeedsRefresh) {
            return;
        }

        if (expNeedsRefresh) {
            expSolver.compute(blockDParams, dt, expD);
            System.arraycopy(blockDParams, 0, cachedBlockDParams, 0, blockDParams.length);
            cachedExpDt = dt;
            expCacheValid = true;
        }

        System.arraycopy(rData, 0, rMatrix.data, 0, d * d);
        System.arraycopy(rtData, 0, rtMatrix.data, 0, d * d);
        multiplyRowMajor(rData, expD.data, d, workMatrix);
        multiply(workMatrix, rtData, d, transitionMatrix);
        System.arraycopy(rData, 0, cachedRData, 0, rData.length);
        System.arraycopy(rtData, 0, cachedRtData, 0, rtData.length);
        cachedBasisDt = dt;
        basisCacheValid = true;
    }

    private void refreshMeanGradientCaches(final double dt) {
        blockParameter.fillBlockDiagonalElements(blockDParams);
        if (!expCacheValid
                || Double.doubleToLongBits(dt) != Double.doubleToLongBits(cachedExpDt)
                || !Arrays.equals(blockDParams, cachedBlockDParams)) {
            expSolver.compute(blockDParams, dt, expD);
            System.arraycopy(blockDParams, 0, cachedBlockDParams, 0, blockDParams.length);
            cachedExpDt = dt;
            expCacheValid = true;
            basisCacheValid = false;
        }
    }

    private static void multiplyRowMajor(final double[] leftRowMajor,
                                         final double[] rightRowMajor,
                                         final int dimension,
                                         final double[][] out) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += leftRowMajor[i * dimension + k] * rightRowMajor[k * dimension + j];
                }
                out[i][j] = sum;
            }
        }
    }

    private static void multiply(final double[][] left,
                                 final double[] rightRowMajor,
                                 final int dimension,
                                 final DenseMatrix64F out) {
        final double[] outData = out.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += left[i][k] * rightRowMajor[k * dimension + j];
                }
                outData[i * dimension + j] = sum;
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

    private void fillTotalUpstreamOnTransition(final double[] stationaryMean,
                                               final double[] dLogL_dF,
                                               final double[] dLogL_df,
                                               final DenseMatrix64F out) {
        final int dimension = getDimension();
        final double[] outData = out.data;
        for (int row = 0; row < dimension; ++row) {
            final int rowOff = row * dimension;
            for (int col = 0; col < dimension; ++col) {
                outData[rowOff + col] =
                        dLogL_dF[rowOff + col] - dLogL_df[row] * stationaryMean[col];
            }
        }
    }

    private void accumulateCompressedGradient(final DenseMatrix64F denseGradient,
                                              final double[] compressedAccumulator) {
        final double[] data = denseGradient.data;
        final int dimension = getDimension();
        final int upperBase = dimension;
        final int lowerBase = dimension + blockParameter.getNum2x2Blocks();
        int blockIndex = 0;

        for (int i = 0; i < dimension; ++i) {
            compressedAccumulator[i] += data[i * dimension + i];
        }
        for (int b = 0; b < blockParameter.getNumBlocks(); ++b) {
            if (blockParameter.getBlockSizes()[b] != 2) {
                continue;
            }
            final int start = blockParameter.getBlockStarts()[b];
            compressedAccumulator[upperBase + blockIndex] += data[start * dimension + (start + 1)];
            compressedAccumulator[lowerBase + blockIndex] += data[(start + 1) * dimension + start];
            blockIndex++;
        }
    }

    private static void addDenseMatrixToArray(final DenseMatrix64F src, final double[][] dest) {
        final int dimension = dest.length;
        final double[] data = src.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                dest[i][j] += data[i * dimension + j];
            }
        }
    }

    private void applyBlockDiagonalTransposeExp(final double[] in,
                                                final double[] out) {
        for (int b = 0; b < blockParameter.getNumBlocks(); ++b) {
            final int start = blockParameter.getBlockStarts()[b];
            final int size = blockParameter.getBlockSizes()[b];
            if (size == 1) {
                out[start] = expD.data[start * getDimension() + start] * in[start];
            } else {
                final int row0 = start * getDimension();
                final int row1 = (start + 1) * getDimension();
                final double e00 = expD.data[row0 + start];
                final double e01 = expD.data[row0 + start + 1];
                final double e10 = expD.data[row1 + start];
                final double e11 = expD.data[row1 + start + 1];
                out[start] = e00 * in[start] + e10 * in[start + 1];
                out[start + 1] = e01 * in[start] + e11 * in[start + 1];
            }
        }
    }

    private void applyBlockDiagonalTransposeExpPrepared(final PreparedBranchBasis prepared,
                                                        final double[] in,
                                                        final double[] out) {
        final double[] expData = prepared.expD.data;
        for (int b = 0; b < blockParameter.getNumBlocks(); ++b) {
            final int start = blockParameter.getBlockStarts()[b];
            final int size = blockParameter.getBlockSizes()[b];
            if (size == 1) {
                out[start] = expData[start * prepared.dimension + start] * in[start];
            } else {
                final int row0 = start * prepared.dimension;
                final int row1 = (start + 1) * prepared.dimension;
                final double e00 = expData[row0 + start];
                final double e01 = expData[row0 + start + 1];
                final double e10 = expData[row1 + start];
                final double e11 = expData[row1 + start + 1];
                out[start] = e00 * in[start] + e10 * in[start + 1];
                out[start + 1] = e01 * in[start] + e11 * in[start + 1];
            }
        }
    }

    private static void fillSymmetricDenseMatrix(final double[][] source, final DenseMatrix64F out) {
        final int dimension = out.numRows;
        final double[] data = out.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                data[i * dimension + j] = 0.5 * (source[i][j] + source[j][i]);
            }
        }
    }

    private static void copyDenseMatrixToArray(final DenseMatrix64F source, final double[][] out) {
        final int dimension = source.numRows;
        final double[] data = source.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                out[i][j] = data[i * dimension + j];
            }
        }
    }

    private static void copyDenseMatrixToFlat(final DenseMatrix64F source, final double[] out) {
        System.arraycopy(source.data, 0, out, 0, source.numRows * source.numCols);
    }

    private static void copySquareMatrixToFlat(final double[][] source, final double[] out) {
        final int dimension = source.length;
        for (int i = 0; i < dimension; ++i) {
            System.arraycopy(source[i], 0, out, i * dimension, dimension);
        }
    }

    private static void fillSymmetricDenseMatrixFlat(final double[] source, final DenseMatrix64F out) {
        final int dimension = out.numRows;
        final double[] data = out.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                data[i * dimension + j] = 0.5 * (source[i * dimension + j] + source[j * dimension + i]);
            }
        }
    }

    private void fillScaledNegativeBlockDMatrix(final double dt,
                                                final double[][] out) {
        final int d = getDimension();
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                out[i][j] = 0.0;
            }
        }
        final int upperOffset = d;
        final int lowerOffset = d + (d - 1);
        for (int i = 0; i < d; ++i) {
            out[i][i] = -dt * blockDParams[i];
            if (i < d - 1) {
                out[i][i + 1] = -dt * blockDParams[upperOffset + i];
                out[i + 1][i] = -dt * blockDParams[lowerOffset + i];
            }
        }
    }

    private static void fillScaledNegativeBlockDMatrix(final double[] blockDParams,
                                                       final double dt,
                                                       final double[][] out,
                                                       final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                out[i][j] = 0.0;
            }
        }
        final int upperOffset = dimension;
        final int lowerOffset = dimension + (dimension - 1);
        for (int i = 0; i < dimension; ++i) {
            out[i][i] = -dt * blockDParams[i];
            if (i < dimension - 1) {
                out[i][i + 1] = -dt * blockDParams[upperOffset + i];
                out[i + 1][i] = -dt * blockDParams[lowerOffset + i];
            }
        }
    }

    private static void multiplyRowMajorVector(final double[] matrix,
                                               final double[] vector,
                                               final double[] out,
                                               final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            double sum = 0.0;
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                sum += matrix[rowOffset + j] * vector[j];
            }
            out[i] = sum;
        }
    }

    private static void fillDenseMatrix(final double[][] source,
                                        final DenseMatrix64F out) {
        final int dimension = out.numRows;
        final double[] data = out.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                data[i * dimension + j] = source[i][j];
            }
        }
    }

    private double copyAndInvertPositiveDefinite(final DenseMatrix64F source,
                                                 final double[][] matrixOut,
                                                 final double[][] inverseOut) {
        final int d = source.numRows;
        final double[] sourceData = source.data;
        double maxAbsDiagonal = 0.0;
        for (int i = 0; i < d; ++i) {
            final double[] matrixRow = matrixOut[i];
            final int iOff = i * d;
            for (int j = 0; j < d; ++j) {
                final double value = 0.5 * (sourceData[iOff + j] + sourceData[j * d + i]);
                matrixRow[j] = value;
                choleskyScratch[iOff + j] = value;
            }
            maxAbsDiagonal = Math.max(maxAbsDiagonal, Math.abs(matrixRow[i]));
        }

        final double pivotFloor = Math.max(
                SPD_JITTER_ABSOLUTE,
                SPD_JITTER_RELATIVE * Math.max(1.0, maxAbsDiagonal));
        int clippedPivotCount = 0;
        double minPivotBeforeFloor = Double.POSITIVE_INFINITY;
        double logDet = 0.0;
        for (int i = 0; i < d; ++i) {
            final int iOff = i * d;
            for (int j = 0; j <= i; ++j) {
                double sum = choleskyScratch[iOff + j];
                for (int k = 0; k < j; ++k) {
                    sum -= choleskyScratch[iOff + k] * choleskyScratch[j * d + k];
                }
                if (i == j) {
                    minPivotBeforeFloor = Math.min(minPivotBeforeFloor, sum);
                    if (sum < pivotFloor) {
                        clippedPivotCount++;
                        sum = pivotFloor;
                    }
                    choleskyScratch[iOff + i] = Math.sqrt(sum);
                    logDet += Math.log(choleskyScratch[iOff + i]);
                } else {
                    final double denominator = Math.max(choleskyScratch[j * d + j], Math.sqrt(pivotFloor));
                    choleskyScratch[iOff + j] = sum / denominator;
                }
            }
        }
        if (Boolean.getBoolean(DEBUG_PD_PIVOT_FLOOR_PROPERTY) && clippedPivotCount > 0) {
            System.err.println(
                    "pdPivotFloorDebug clipped=" + clippedPivotCount
                            + " dim=" + d
                            + " minPivotBeforeFloor=" + minPivotBeforeFloor
                            + " pivotFloor=" + pivotFloor
                            + " maxAbsDiagonal=" + maxAbsDiagonal);
        }

        for (int col = 0; col < d; ++col) {
            for (int row = 0; row < d; ++row) {
                double sum = (row == col) ? 1.0 : 0.0;
                final int rowOff = row * d;
                for (int k = 0; k < row; ++k) {
                    sum -= choleskyScratch[rowOff + k] * lowerInverseScratch[k * d + col];
                }
                final double denominator = Math.max(choleskyScratch[rowOff + row], Math.sqrt(pivotFloor));
                lowerInverseScratch[rowOff + col] = sum / denominator;
            }
        }

        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += lowerInverseScratch[k * d + i] * lowerInverseScratch[k * d + j];
                }
                inverseOut[i][j] = sum;
            }
        }
        return 2.0 * logDet;
    }

    private double copyAndInvertPositiveDefiniteFlat(final DenseMatrix64F source,
                                                     final double[] matrixOut,
                                                     final double[] inverseOut) {
        return copyAndInvertPositiveDefiniteFlat(
                source, matrixOut, inverseOut, choleskyScratch, lowerInverseScratch);
    }

    private static double copyAndInvertPositiveDefiniteFlat(final DenseMatrix64F source,
                                                            final double[] matrixOut,
                                                            final double[] inverseOut,
                                                            final double[] choleskyScratch,
                                                            final double[] lowerInverseScratch) {
        final int d = source.numRows;
        final double[] sourceData = source.data;
        double maxAbsDiagonal = 0.0;
        for (int i = 0; i < d; ++i) {
            final int iOff = i * d;
            for (int j = 0; j < d; ++j) {
                final double value = 0.5 * (sourceData[iOff + j] + sourceData[j * d + i]);
                matrixOut[iOff + j] = value;
                choleskyScratch[iOff + j] = value;
            }
            maxAbsDiagonal = Math.max(maxAbsDiagonal, Math.abs(matrixOut[iOff + i]));
        }

        final double pivotFloor = Math.max(
                SPD_JITTER_ABSOLUTE,
                SPD_JITTER_RELATIVE * Math.max(1.0, maxAbsDiagonal));
        int clippedPivotCount = 0;
        double minPivotBeforeFloor = Double.POSITIVE_INFINITY;
        double logDet = 0.0;
        for (int i = 0; i < d; ++i) {
            final int iOff = i * d;
            for (int j = 0; j <= i; ++j) {
                double sum = choleskyScratch[iOff + j];
                for (int k = 0; k < j; ++k) {
                    sum -= choleskyScratch[iOff + k] * choleskyScratch[j * d + k];
                }
                if (i == j) {
                    minPivotBeforeFloor = Math.min(minPivotBeforeFloor, sum);
                    if (sum < pivotFloor) {
                        clippedPivotCount++;
                        sum = pivotFloor;
                    }
                    choleskyScratch[iOff + i] = Math.sqrt(sum);
                    logDet += Math.log(choleskyScratch[iOff + i]);
                } else {
                    final double denominator = Math.max(choleskyScratch[j * d + j], Math.sqrt(pivotFloor));
                    choleskyScratch[iOff + j] = sum / denominator;
                }
            }
        }
        if (Boolean.getBoolean(DEBUG_PD_PIVOT_FLOOR_PROPERTY) && clippedPivotCount > 0) {
            System.err.println(
                    "pdPivotFloorDebug clipped=" + clippedPivotCount
                            + " dim=" + d
                            + " minPivotBeforeFloor=" + minPivotBeforeFloor
                            + " pivotFloor=" + pivotFloor
                            + " maxAbsDiagonal=" + maxAbsDiagonal);
        }

        for (int col = 0; col < d; ++col) {
            for (int row = 0; row < d; ++row) {
                double sum = (row == col) ? 1.0 : 0.0;
                final int rowOff = row * d;
                for (int k = 0; k < row; ++k) {
                    sum -= choleskyScratch[rowOff + k] * lowerInverseScratch[k * d + col];
                }
                final double denominator = Math.max(choleskyScratch[rowOff + row], Math.sqrt(pivotFloor));
                lowerInverseScratch[rowOff + col] = sum / denominator;
            }
        }

        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += lowerInverseScratch[k * d + i] * lowerInverseScratch[k * d + j];
                }
                inverseOut[i * d + j] = sum;
            }
        }
        return 2.0 * logDet;
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

    private static void transposeInPlace(final DenseMatrix64F matrix) {
        final int dimension = matrix.numRows;
        final double[] data = matrix.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = i + 1; j < dimension; ++j) {
                final int ij = i * dimension + j;
                final int ji = j * dimension + i;
                final double tmp = data[ij];
                data[ij] = data[ji];
                data[ji] = tmp;
            }
        }
    }

    private static boolean isFinite(final double[] values) {
        for (double value : values) {
            if (!Double.isFinite(value)) {
                return false;
            }
        }
        return true;
    }

    private static boolean isFinite(final double[][] values) {
        for (double[] row : values) {
            if (!isFinite(row)) {
                return false;
            }
        }
        return true;
    }

}
