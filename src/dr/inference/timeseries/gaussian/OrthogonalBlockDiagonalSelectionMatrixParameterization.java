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

/**
 * Orthogonal block-diagonal parametrization using block-space helpers for the
 * forward path and the orthogonal native backward path.
 */
public final class OrthogonalBlockDiagonalSelectionMatrixParameterization extends DenseSelectionMatrixParameterization {

    private static final double SPD_JITTER_RELATIVE = 1.0e-14;
    private static final double SPD_JITTER_ABSOLUTE = 1.0e-14;
    private static final String TRANSPOSE_NATIVE_FRECHET_OUTPUT_PROPERTY =
            "beast.experimental.transposeNativeFrechetOutput";
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
    private final DenseMatrix64F expD;
    private final DenseMatrix64F rMatrix;
    private final DenseMatrix64F rtMatrix;
    private final DenseMatrix64F transitionMatrix;
    private final DenseMatrix64F qMatrix;
    private final DenseMatrix64F qDBasis;
    private final DenseMatrix64F stationaryCovDBasis;
    private final DenseMatrix64F transitionCovDBasis;
    private final DenseMatrix64F transitionCovariance;
    private final double[][] transitionMatrixArrayScratch;
    private final double[][] transitionCovarianceArrayScratch;
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
    private final double[][] choleskyScratch;
    private final double[][] lowerInverseScratch;
    private final CanonicalTransitionAdjointUtils.Workspace canonicalAdjointWorkspace;
    private final double[] tempVector1;
    private final double[] tempVector2;

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
        this.expD = new DenseMatrix64F(d, d);
        this.rMatrix = new DenseMatrix64F(d, d);
        this.rtMatrix = new DenseMatrix64F(d, d);
        this.transitionMatrix = new DenseMatrix64F(d, d);
        this.qMatrix = new DenseMatrix64F(d, d);
        this.qDBasis = new DenseMatrix64F(d, d);
        this.stationaryCovDBasis = new DenseMatrix64F(d, d);
        this.transitionCovDBasis = new DenseMatrix64F(d, d);
        this.transitionCovariance = new DenseMatrix64F(d, d);
        this.transitionMatrixArrayScratch = new double[d][d];
        this.transitionCovarianceArrayScratch = new double[d][d];
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
        this.choleskyScratch = new double[d][d];
        this.lowerInverseScratch = new double[d][d];
        this.canonicalAdjointWorkspace = new CanonicalTransitionAdjointUtils.Workspace(d);
        this.tempVector1 = new double[d];
        this.tempVector2 = new double[d];
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

    public void fillCanonicalLocalAdjoints(final MatrixParameterInterface diffusionMatrix,
                                           final double[] stationaryMean,
                                           final double dt,
                                           final CanonicalBranchMessageContribution contribution,
                                           final CanonicalLocalTransitionAdjoints out) {
        refreshBasisCaches(dt);
        fillTransitionCovarianceMatrix(diffusionMatrix, dt, transitionCovariance);
        copyDenseMatrixToArray(transitionMatrix, transitionMatrixArrayScratch);
        fillTransitionOffset(stationaryMean, transitionOffsetScratch);

        final double[][] precision = workMatrix;
        copyAndInvertPositiveDefinite(transitionCovariance, transitionCovarianceArrayScratch, precision);
        CanonicalTransitionAdjointUtils.fillFromMoments(
                precision,
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

        final double logDet = copyAndInvertPositiveDefinite(
                transitionCovariance,
                transitionCovarianceArrayScratch,
                out.precisionYY);

        // J_yx = -P F
        for (int i = 0; i < d; ++i) {
            final double[] pRow = out.precisionYY[i];
            final double[] jyxRow = out.precisionYX[i];
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += pRow[k] * transitionData[k * d + j];
                }
                jyxRow[j] = -sum;
            }
        }

        // J_xy = J_yx^T
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                out.precisionXY[i][j] = out.precisionYX[j][i];
            }
        }

        // J_xx = F^T P F = -(F^T J_yx)
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum -= transitionData[k * d + i] * out.precisionYX[k][j];
                }
                out.precisionXX[i][j] = sum;
            }
        }

        // h_y = P f
        for (int i = 0; i < d; ++i) {
            double sum = 0.0;
            final double[] pRow = out.precisionYY[i];
            for (int j = 0; j < d; ++j) {
                sum += pRow[j] * transitionOffset[j];
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
        refreshBasisCaches(dt);
        accumulateNativeGradientFromTransitionCached(
                dt, stationaryMean, dLogL_dF, dLogL_df, compressedDAccumulator, rotationAccumulator);
    }

    public void accumulateNativeGradientFromCovarianceStationary(final MatrixParameterInterface diffusionMatrix,
                                                                 final double dt,
                                                                 final double[][] dLogL_dV,
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
        copyDenseMatrixToArray(transitionMatrix, transitionMatrixArrayScratch);
        fillTransitionOffset(stationaryMean, transitionOffsetScratch);

        final double[][] precision = workMatrix;
        copyAndInvertPositiveDefinite(transitionCovariance, transitionCovarianceArrayScratch, precision);
        CanonicalTransitionAdjointUtils.fillFromMoments(
                precision,
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
                                                              final double[][] dLogL_dF,
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
            if (Boolean.getBoolean(TRANSPOSE_NATIVE_FRECHET_OUTPUT_PROPERTY)) {
                transposeInPlace(gradD);
            }
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
            if (Boolean.getBoolean(TRANSPOSE_NATIVE_FRECHET_OUTPUT_PROPERTY)) {
                transposeInPlace(gradD);
            }
            if (!isFinite(gradD.data)) {
                fillScaledNegativeBlockDMatrix(dt, scaledNegativeBlockDScratch);
                copyDenseMatrixToArray(upstreamFD, denseAdjointScratch);
                MatrixExponentialUtils.adjointExp(
                        scaledNegativeBlockDScratch,
                        denseAdjointScratch,
                        denseAdjointScratch);
                fillDenseMatrix(denseAdjointScratch, gradD);
                if (Boolean.getBoolean(TRANSPOSE_NATIVE_FRECHET_OUTPUT_PROPERTY)) {
                    transposeInPlace(gradD);
                }
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

    private void accumulateNativeGradientFromCovarianceStationaryCached(final double dt,
                                                                        final double[][] dLogL_dV,
                                                                        final double[] compressedDAccumulator,
                                                                        final double[][] rotationAccumulator) {
        final boolean forceDenseAdjointExp = Boolean.getBoolean(NATIVE_FORCE_DENSE_ADJOINT_EXP_PROPERTY);
        fillSymmetricDenseMatrix(dLogL_dV, gV);

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
            if (Boolean.getBoolean(TRANSPOSE_NATIVE_FRECHET_OUTPUT_PROPERTY)) {
                transposeInPlace(gradD);
            }
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
            if (Boolean.getBoolean(TRANSPOSE_NATIVE_FRECHET_OUTPUT_PROPERTY)) {
                transposeInPlace(gradD);
            }
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

    private void refreshBasisCaches(final double dt) {
        final int d = getDimension();
        orthogonalRotation.fillOrthogonalMatrix(rData);
        orthogonalRotation.fillOrthogonalTranspose(rtData);
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
        System.arraycopy(rData, 0, rMatrix.data, 0, d * d);
        System.arraycopy(rtData, 0, rtMatrix.data, 0, d * d);
        blockParameter.fillBlockDiagonalElements(blockDParams);
        expSolver.compute(blockDParams, dt, expD);
        multiplyRowMajor(rData, expD.data, d, workMatrix);
        multiply(workMatrix, rtData, d, transitionMatrix);
    }

    private void refreshMeanGradientCaches(final double dt) {
        blockParameter.fillBlockDiagonalElements(blockDParams);
        expSolver.compute(blockDParams, dt, expD);
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
                                               final double[][] dLogL_dF,
                                               final double[] dLogL_df,
                                               final DenseMatrix64F out) {
        final int dimension = getDimension();
        final double[] outData = out.data;
        for (int row = 0; row < dimension; ++row) {
            for (int col = 0; col < dimension; ++col) {
                outData[row * dimension + col] =
                        dLogL_dF[row][col] - dLogL_df[row] * stationaryMean[col];
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
            final double[] cholRow = choleskyScratch[i];
            final int rowOffset = i * d;
            for (int j = 0; j < d; ++j) {
                final double value = 0.5 * (sourceData[rowOffset + j] + sourceData[j * d + i]);
                matrixRow[j] = value;
                cholRow[j] = value;
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
            final double[] cholRow = choleskyScratch[i];
            for (int j = 0; j <= i; ++j) {
                double sum = cholRow[j];
                for (int k = 0; k < j; ++k) {
                    sum -= cholRow[k] * choleskyScratch[j][k];
                }
                if (i == j) {
                    minPivotBeforeFloor = Math.min(minPivotBeforeFloor, sum);
                    if (sum < pivotFloor) {
                        clippedPivotCount++;
                        sum = pivotFloor;
                    }
                    cholRow[j] = Math.sqrt(sum);
                    logDet += Math.log(cholRow[j]);
                } else {
                    final double denominator = Math.max(choleskyScratch[j][j], Math.sqrt(pivotFloor));
                    cholRow[j] = sum / denominator;
                }
            }
            for (int j = i + 1; j < d; ++j) {
                cholRow[j] = 0.0;
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
                for (int k = 0; k < row; ++k) {
                    sum -= choleskyScratch[row][k] * lowerInverseScratch[k][col];
                }
                final double denominator = Math.max(choleskyScratch[row][row], Math.sqrt(pivotFloor));
                lowerInverseScratch[row][col] = sum / denominator;
            }
        }

        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += lowerInverseScratch[k][i] * lowerInverseScratch[k][j];
                }
                inverseOut[i][j] = sum;
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
