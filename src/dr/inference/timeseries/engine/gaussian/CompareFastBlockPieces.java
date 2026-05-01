package dr.inference.timeseries.engine.gaussian;

import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.BlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.core.UniformTimeGrid;
import dr.inference.timeseries.gaussian.BlockDiagonalSelectionMatrixParameterization;
import dr.inference.timeseries.gaussian.GaussianObservationModel;
import dr.inference.timeseries.gaussian.OUTimeSeriesProcessAdapter;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;

public class CompareFastBlockPieces {

    private static MatrixParameter makeMatrix(String name, double[][] values) {
        MatrixParameter mp = new MatrixParameter(name, values.length, values[0].length);
        for (int i = 0; i < values.length; i++) {
            for (int j = 0; j < values[0].length; j++) {
                mp.setParameterValue(i, j, values[i][j]);
            }
        }
        return mp;
    }

    private static final class Model {
        final OUProcessModel process;
        final AnalyticalKalmanGradientEngine analytical;
        final KalmanSmootherEngine smoother;

        private Model(OUProcessModel process,
                      AnalyticalKalmanGradientEngine analytical,
                      KalmanSmootherEngine smoother) {
            this.process = process;
            this.analytical = analytical;
            this.smoother = smoother;
        }
    }

    private static Model makeModel() {
        final double[][] rotationValues = {{1.0, 0.25}, {-0.15, 0.9}};
        final double[][] diffusion = {{1.2, 0.1}, {0.1, 1.6}};
        final double[][] initCov = {{1.0, 0.2}, {0.2, 1.3}};
        final double[][] noise = {{0.6, 0.05}, {0.05, 0.8}};
        final double[] y1 = {0.7, -0.4, 1.1, 0.2};
        final double[] y2 = {-0.3, 0.8, -1.2, 0.5};
        final double dt = 0.15;

        MatrixParameter rotation = makeMatrix("R", rotationValues);
        Parameter scalarBlock = new Parameter.Default(0);
        Parameter rhoParam = new Parameter.Default(0.02);
        Parameter thetaParam = new Parameter.Default(0.0);
        Parameter tParam = new Parameter.Default(0.0);
        MatrixParameter drift = new BlockDiagonalPolarStableMatrixParameter(
                "Ablock", rotation, scalarBlock, rhoParam, thetaParam, tParam);
        MatrixParameter diffusionParam = makeMatrix("Q", diffusion);
        Parameter mean = new Parameter.Default(0.25);
        MatrixParameter initCovParam = makeMatrix("P0", initCov);
        OUProcessModel process = new OUProcessModel("ou2dblock", 2, drift, diffusionParam, mean, initCovParam,
                OUProcessModel.CovarianceGradientMethod.STATIONARY_LYAPUNOV);

        MatrixParameter H = makeMatrix("H", new double[][]{{1, 0}, {0, 1}});
        MatrixParameter R = makeMatrix("Robs", noise);
        MatrixParameter Y = makeMatrix("Y", new double[][]{y1, y2});
        GaussianObservationModel obs = new GaussianObservationModel("obs", 2, H, R, Y);
        TimeGrid grid = new UniformTimeGrid(y1.length, 0.0, dt);
        GaussianTransitionRepresentation rep = new OUTimeSeriesProcessAdapter(process)
                .getRepresentation(GaussianTransitionRepresentation.class);
        KalmanSmootherEngine smoother = new KalmanSmootherEngine(rep, obs, grid);
        AnalyticalKalmanGradientEngine analytical = new AnalyticalKalmanGradientEngine(
                smoother, new SelectionMatrixGradientFormula(process.getDriftMatrix(), 2));
        return new Model(process, analytical, smoother);
    }

    public static void main(String[] args) {
        final Model model = makeModel();
        final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParam =
                (AbstractBlockDiagonalTwoByTwoMatrixParameter) model.process.getDriftMatrix();
        final BlockDiagonalSelectionMatrixParameterization parameterization =
                (BlockDiagonalSelectionMatrixParameterization) model.process.getSelectionMatrixParameterization();

        final BranchSmootherStats[] smootherStats = model.smoother.getSmootherStats();
        final ForwardTrajectory trajectory = model.smoother.getTrajectory();
        final int T = trajectory.timeCount;
        final int d = 2;

        final double[] mu = new double[d];
        model.process.getInitialMean(mu);

        final double[] fastTransitionD = new double[blockParam.getCompressedDDimension()];
        final double[] fastCovarianceD = new double[blockParam.getCompressedDDimension()];
        final double[][] fastTransitionR = new double[d][d];
        final double[][] fastCovarianceR = new double[d][d];

        final double[] denseTransitionA = new double[d * d];
        final double[] denseCovarianceA = new double[d * d];

        final double[] meanResidual = new double[d];
        final double[][] crossCov = new double[d][d];
        final double[][] branchMat = new double[d][d];
        final double[][] stepCovInv = new double[d][d];
        final double[][] temp = new double[d][d];
        final double[][] temp2 = new double[d][d];
        final double[][] dLogL_dF = new double[d][d];
        final double[] dLogL_df = new double[d];
        final double[][] residualSecondMoment = new double[d][d];
        final double[][] dLogL_dV = new double[d][d];

        for (int t = 0; t < T - 1; ++t) {
            final BranchSmootherStats curr = smootherStats[t];
            final BranchSmootherStats next = smootherStats[t + 1];
            final double[][] F_t = trajectory.transitionMatrices[t];
            final double[] f_t = trajectory.transitionOffsets[t];
            final double dt = model.smoother.getTimeGrid().getDelta(t, t + 1);

            KalmanLikelihoodEngine.multiplyMatrixMatrixTransposedRight(
                    next.smoothedCovariance, curr.smootherGain, crossCov);
            KalmanLikelihoodEngine.multiplyMatrixVector(F_t, curr.smoothedMean, meanResidual);
            for (int i = 0; i < d; ++i) {
                meanResidual[i] = next.smoothedMean[i] - meanResidual[i] - f_t[i];
            }

            KalmanLikelihoodEngine.copyMatrix(crossCov, branchMat);
            KalmanLikelihoodEngine.multiplyMatrixMatrix(F_t, curr.smoothedCovariance, temp);
            for (int i = 0; i < d; ++i) {
                for (int j = 0; j < d; ++j) {
                    branchMat[i][j] -= temp[i][j];
                    branchMat[i][j] += meanResidual[i] * curr.smoothedMean[j];
                }
            }

            KalmanLikelihoodEngine.copyMatrix(trajectory.stepCovariances[t], stepCovInv);
            final KalmanLikelihoodEngine.CholeskyFactor stepChol =
                    KalmanLikelihoodEngine.cholesky(stepCovInv);
            KalmanLikelihoodEngine.invertPositiveDefiniteFromCholesky(stepCovInv, stepChol);

            KalmanLikelihoodEngine.multiplyMatrixMatrix(stepCovInv, branchMat, dLogL_dF);
            KalmanLikelihoodEngine.multiplyMatrixVector(stepCovInv, meanResidual, dLogL_df);

            parameterization.accumulateNativeGradientFromTransition(
                    dt, mu, dLogL_dF, dLogL_df, fastTransitionD, fastTransitionR);
            model.process.accumulateSelectionGradient(dt, dLogL_dF, dLogL_df, denseTransitionA);

            KalmanLikelihoodEngine.copyMatrix(next.smoothedCovariance, residualSecondMoment);
            KalmanLikelihoodEngine.multiplyMatrixMatrixTransposedRight(F_t, crossCov, temp);
            for (int i = 0; i < d; ++i) {
                for (int j = 0; j < d; ++j) {
                    residualSecondMoment[i][j] -= temp[i][j];
                }
            }
            KalmanLikelihoodEngine.multiplyMatrixMatrixTransposedRight(crossCov, F_t, temp);
            for (int i = 0; i < d; ++i) {
                for (int j = 0; j < d; ++j) {
                    residualSecondMoment[i][j] -= temp[i][j];
                }
            }
            KalmanLikelihoodEngine.multiplyMatrixMatrix(F_t, curr.smoothedCovariance, temp);
            KalmanLikelihoodEngine.multiplyMatrixMatrixTransposedRight(temp, F_t, temp2);
            for (int i = 0; i < d; ++i) {
                for (int j = 0; j < d; ++j) {
                    residualSecondMoment[i][j] += temp2[i][j];
                    residualSecondMoment[i][j] += meanResidual[i] * meanResidual[j];
                }
            }

            KalmanLikelihoodEngine.multiplyMatrixMatrix(stepCovInv, residualSecondMoment, temp);
            KalmanLikelihoodEngine.multiplyMatrixMatrix(temp, stepCovInv, dLogL_dV);
            for (int i = 0; i < d; ++i) {
                for (int j = 0; j < d; ++j) {
                    dLogL_dV[i][j] = 0.5 * (dLogL_dV[i][j] - stepCovInv[i][j]);
                }
            }

            parameterization.accumulateNativeGradientFromCovarianceStationary(
                    model.process.getDiffusionMatrix(), dt, dLogL_dV, fastCovarianceD, fastCovarianceR);
            model.process.accumulateSelectionGradientFromCovariance(dt, dLogL_dV, denseCovarianceA);
        }

        final double[] denseTransitionNative = pullBack(blockParam, denseTransitionA);
        final double[] denseCovarianceNative = pullBack(blockParam, denseCovarianceA);
        final double[] fastTransitionNative = project(blockParam, fastTransitionD, fastTransitionR);
        final double[] fastCovarianceNative = project(blockParam, fastCovarianceD, fastCovarianceR);

        printComparison("transition", fastTransitionNative, denseTransitionNative);
        printComparison("covariance", fastCovarianceNative, denseCovarianceNative);
    }

    private static double[] project(AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                    double[] compressedBlockGradient,
                                    double[][] rotationGradient) {
        final double[] nativeGradient = new double[blockParameter.getBlockDiagonalNParameters()];
        blockParameter.chainGradient(compressedBlockGradient, nativeGradient);
        final double[] out = new double[nativeGradient.length + rotationGradient.length * rotationGradient.length];
        System.arraycopy(nativeGradient, 0, out, 0, nativeGradient.length);
        flattenColumnMajor(rotationGradient, out, nativeGradient.length);
        return out;
    }

    private static double[] pullBack(AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                     double[] denseGradient) {
        final int d = blockParameter.getRowDimension();
        final double[][] gradientA = new double[d][d];
        final double[] rData = new double[d * d];
        final double[] rinvData = new double[d * d];
        final double[] rawBlockSource = new double[blockParameter.getCompressedDDimension()];
        final double[] nativeGradient = new double[blockParameter.getBlockDiagonalNParameters()];
        final double[][] gradientD = new double[d][d];
        final double[][] gradientR = new double[d][d];
        final double[][] aMatrix = new double[d][d];

        fillMatrixFromColumnMajorGradient(gradientA, denseGradient, d);
        blockParameter.fillRAndRinv(rData, rinvData);
        copyMatrixParameter(blockParameter, aMatrix);
        computeGradientWrtBlockDiagonalBasis(gradientA, rData, rinvData, gradientD);
        computeGradientWrtRotationMatrix(gradientA, aMatrix, rinvData, gradientR);
        compressActiveBlockGradient(blockParameter, gradientD, rawBlockSource);
        blockParameter.chainGradient(rawBlockSource, nativeGradient);

        final double[] out = new double[nativeGradient.length + d * d];
        System.arraycopy(nativeGradient, 0, out, 0, nativeGradient.length);
        flattenColumnMajor(gradientR, out, nativeGradient.length);
        return out;
    }

    private static void printComparison(String label, double[] fast, double[] dense) {
        System.out.println(label);
        for (int i = 0; i < fast.length; ++i) {
            System.out.println("  idx=" + i + " fast=" + fast[i] + " dense=" + dense[i] + " diff=" + (fast[i] - dense[i]));
        }
    }

    private static void fillMatrixFromColumnMajorGradient(final double[][] out,
                                                          final double[] gradient,
                                                          final int d) {
        for (int col = 0; col < d; ++col) {
            for (int row = 0; row < d; ++row) {
                out[row][col] = gradient[col * d + row];
            }
        }
    }

    private static void computeGradientWrtBlockDiagonalBasis(final double[][] gradientA,
                                                             final double[] rData,
                                                             final double[] rinvData,
                                                             final double[][] out) {
        final int d = gradientA.length;
        final double[][] temp = new double[d][d];
        final double[][] rinvTranspose = new double[d][d];

        for (int row = 0; row < d; ++row) {
            for (int col = 0; col < d; ++col) {
                rinvTranspose[row][col] = rinvData[col * d + row];
            }
        }

        multiplyTransposeLeft(rData, d, gradientA, temp);
        multiply(temp, rinvTranspose, out);
    }

    private static void computeGradientWrtRotationMatrix(final double[][] gradientA,
                                                         final double[][] aMatrix,
                                                         final double[] rinvData,
                                                         final double[][] out) {
        final int d = gradientA.length;
        final double[][] aTranspose = new double[d][d];
        final double[][] temp = new double[d][d];
        final double[][] rinvTranspose = new double[d][d];

        transpose(aMatrix, aTranspose);
        multiply(gradientA, aTranspose, temp);
        multiply(aTranspose, gradientA, out);

        for (int row = 0; row < d; ++row) {
            for (int col = 0; col < d; ++col) {
                rinvTranspose[row][col] = rinvData[col * d + row];
            }
        }

        subtractInPlace(temp, out);
        multiply(temp, rinvTranspose, out);
    }

    private static void compressActiveBlockGradient(final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                                    final double[][] gradientD,
                                                    final double[] out) {
        final int d = gradientD.length;
        int blockIndex = 0;
        final int upperBase = d;
        final int lowerBase = d + blockParameter.getNum2x2Blocks();

        for (int i = 0; i < d; ++i) {
            out[i] = gradientD[i][i];
        }

        for (int b = 0; b < blockParameter.getNumBlocks(); ++b) {
            if (blockParameter.getBlockSizes()[b] != 2) {
                continue;
            }
            final int start = blockParameter.getBlockStarts()[b];
            out[upperBase + blockIndex] = gradientD[start][start + 1];
            out[lowerBase + blockIndex] = gradientD[start + 1][start];
            blockIndex++;
        }
    }

    private static void copyMatrixParameter(final MatrixParameter parameter, final double[][] out) {
        final int d = out.length;
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                out[i][j] = parameter.getParameterValue(i, j);
            }
        }
    }

    private static void flattenColumnMajor(final double[][] matrix, final double[] out, final int offset) {
        final int d = matrix.length;
        int index = offset;
        for (int col = 0; col < d; ++col) {
            for (int row = 0; row < d; ++row) {
                out[index++] = matrix[row][col];
            }
        }
    }

    private static void transpose(final double[][] in, final double[][] out) {
        final int d = in.length;
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                out[i][j] = in[j][i];
            }
        }
    }

    private static void multiplyTransposeLeft(final double[] leftRowMajor,
                                              final int d,
                                              final double[][] right,
                                              final double[][] out) {
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += leftRowMajor[k * d + i] * right[k][j];
                }
                out[i][j] = sum;
            }
        }
    }

    private static void multiply(final double[][] left,
                                 final double[][] right,
                                 final double[][] out) {
        final int d = left.length;
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += left[i][k] * right[k][j];
                }
                out[i][j] = sum;
            }
        }
    }

    private static void subtractInPlace(final double[][] leftMinusRightTarget,
                                        final double[][] right) {
        final int d = leftMinusRightTarget.length;
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                leftMinusRightTarget[i][j] -= right[i][j];
            }
        }
    }
}
