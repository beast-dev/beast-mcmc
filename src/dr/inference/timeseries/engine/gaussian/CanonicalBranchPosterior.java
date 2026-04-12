package dr.inference.timeseries.engine.gaussian;

import dr.inference.timeseries.representation.CanonicalGaussianState;
import dr.inference.timeseries.representation.CanonicalGaussianUtils;

/**
 * Joint posterior statistics for a single branch pair {@code (x_t, x_{t+1})}.
 *
 * <p>The canonical form is the primary representation. Moment-form accessors are
 * derived from it on demand for branch-gradient calculations.
 */
final class CanonicalBranchPosterior {

    private final int stateDimension;
    private final CanonicalGaussianState pairState;
    private final double[] jointMean;
    private final double[][] jointCovariance;
    private final double[] currentMeanScratch;
    private final double[] nextMeanScratch;
    private final double[][] currentCovarianceScratch;
    private final double[][] crossCovarianceScratch;
    private final double[][] nextCovarianceScratch;
    private final double[][] currentSecondMomentScratch;
    private final double[][] crossSecondMomentScratch;
    private final double[][] nextSecondMomentScratch;

    CanonicalBranchPosterior(final int stateDimension) {
        this.stateDimension = stateDimension;
        this.pairState = new CanonicalGaussianState(2 * stateDimension);
        this.jointMean = new double[2 * stateDimension];
        this.jointCovariance = new double[2 * stateDimension][2 * stateDimension];
        this.currentMeanScratch = new double[stateDimension];
        this.nextMeanScratch = new double[stateDimension];
        this.currentCovarianceScratch = new double[stateDimension][stateDimension];
        this.crossCovarianceScratch = new double[stateDimension][stateDimension];
        this.nextCovarianceScratch = new double[stateDimension][stateDimension];
        this.currentSecondMomentScratch = new double[stateDimension][stateDimension];
        this.crossSecondMomentScratch = new double[stateDimension][stateDimension];
        this.nextSecondMomentScratch = new double[stateDimension][stateDimension];
    }

    void fillFromMoments(final double[] currentMean,
                         final double[][] currentCovariance,
                         final double[] nextMean,
                         final double[][] nextCovariance,
                         final double[][] crossCovariance) {
        final int d = stateDimension;
        final double[] jointMeanLocal = jointMean;
        final double[][] jointCovarianceLocal = jointCovariance;

        for (int i = 0; i < d; ++i) {
            jointMeanLocal[i] = currentMean[i];
            jointMeanLocal[d + i] = nextMean[i];
        }
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                jointCovarianceLocal[i][j] = currentCovariance[i][j];
                jointCovarianceLocal[i][d + j] = crossCovariance[j][i];
                jointCovarianceLocal[d + i][j] = crossCovariance[i][j];
                jointCovarianceLocal[d + i][d + j] = nextCovariance[i][j];
            }
        }

        CanonicalGaussianUtils.fillStateFromMoments(jointMeanLocal, jointCovarianceLocal, pairState);
    }

    CanonicalGaussianState getPairState() {
        return pairState;
    }

    void fillFromCanonicalPairState(final CanonicalGaussianState source) {
        for (int i = 0; i < pairState.getDimension(); ++i) {
            pairState.information[i] = source.information[i];
            for (int j = 0; j < pairState.getDimension(); ++j) {
                pairState.precision[i][j] = source.precision[i][j];
            }
        }
        pairState.logNormalizer = source.logNormalizer;
    }

    void extractMoments(final double[] currentMeanOut,
                        final double[][] currentCovarianceOut,
                        final double[] nextMeanOut,
                        final double[][] nextCovarianceOut,
                        final double[][] crossCovarianceOut) {
        CanonicalGaussianUtils.fillMomentsFromCanonical(pairState, jointMean, jointCovariance);
        final int d = stateDimension;
        for (int i = 0; i < d; ++i) {
            currentMeanOut[i] = jointMean[i];
            nextMeanOut[i] = jointMean[d + i];
            for (int j = 0; j < d; ++j) {
                currentCovarianceOut[i][j] = jointCovariance[i][j];
                crossCovarianceOut[i][j] = jointCovariance[d + i][j];
                nextCovarianceOut[i][j] = jointCovariance[d + i][d + j];
            }
        }
    }

    void fillSufficientStatistics(final double[] currentMeanOut,
                                  final double[] nextMeanOut,
                                  final double[][] currentSecondMomentOut,
                                  final double[][] crossSecondMomentOut,
                                  final double[][] nextSecondMomentOut) {
        CanonicalGaussianUtils.fillMomentsFromCanonical(pairState, jointMean, jointCovariance);
        final int d = stateDimension;
        for (int i = 0; i < d; ++i) {
            final double currentMeanI = jointMean[i];
            final double nextMeanI = jointMean[d + i];
            currentMeanOut[i] = currentMeanI;
            nextMeanOut[i] = nextMeanI;
            for (int j = 0; j < d; ++j) {
                final double currentMeanJ = jointMean[j];
                final double nextMeanJ = jointMean[d + j];
                currentSecondMomentOut[i][j] = jointCovariance[i][j] + currentMeanI * currentMeanJ;
                crossSecondMomentOut[i][j] = jointCovariance[d + i][j] + nextMeanI * currentMeanJ;
                nextSecondMomentOut[i][j] = jointCovariance[d + i][d + j] + nextMeanI * nextMeanJ;
            }
        }
    }

    void fillTransitionSufficientStatistics(final double[][] transitionMatrix,
                                            final double[] transitionOffset,
                                            final double[] meanResidualOut,
                                            final double[][] residualCrossMomentOut,
                                            final double[][] residualSecondMomentOut) {
        extractMoments(
                currentMeanScratch,
                currentCovarianceScratch,
                nextMeanScratch,
                nextCovarianceScratch,
                crossCovarianceScratch);
        final int d = stateDimension;
        for (int i = 0; i < d; ++i) {
            double predictedMean = transitionOffset[i];
            for (int k = 0; k < d; ++k) {
                predictedMean += transitionMatrix[i][k] * currentMeanScratch[k];
            }
            meanResidualOut[i] = nextMeanScratch[i] - predictedMean;
        }

        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                residualCrossMomentOut[i][j] = crossCovarianceScratch[i][j];
                residualSecondMomentOut[i][j] = nextCovarianceScratch[i][j];
            }
        }
        multiplyMatrixMatrix(transitionMatrix, currentCovarianceScratch, currentSecondMomentScratch);
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                residualCrossMomentOut[i][j] -= currentSecondMomentScratch[i][j];
                residualCrossMomentOut[i][j] += meanResidualOut[i] * currentMeanScratch[j];
            }
        }
        multiplyMatrixMatrixTransposedRight(transitionMatrix, crossCovarianceScratch, currentSecondMomentScratch);
        subtractInPlace(residualSecondMomentOut, currentSecondMomentScratch);
        multiplyMatrixMatrixTransposedRight(crossCovarianceScratch, transitionMatrix, currentSecondMomentScratch);
        subtractInPlace(residualSecondMomentOut, currentSecondMomentScratch);
        multiplyMatrixMatrix(transitionMatrix, currentCovarianceScratch, currentSecondMomentScratch);
        multiplyMatrixMatrixTransposedRight(currentSecondMomentScratch, transitionMatrix, nextSecondMomentScratch);
        addInPlace(residualSecondMomentOut, nextSecondMomentScratch);
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                residualSecondMomentOut[i][j] += meanResidualOut[i] * meanResidualOut[j];
            }
        }
    }

    void fillLocalMessageContribution(final CanonicalBranchMessageContribution out) {
        fillSufficientStatistics(
                currentMeanScratch,
                nextMeanScratch,
                currentSecondMomentScratch,
                crossSecondMomentScratch,
                nextSecondMomentScratch);
        final int d = stateDimension;
        for (int i = 0; i < d; ++i) {
            out.dLogL_dInformationX[i] = currentMeanScratch[i];
            out.dLogL_dInformationY[i] = nextMeanScratch[i];
            for (int j = 0; j < d; ++j) {
                out.dLogL_dPrecisionXX[i][j] = -0.5 * currentSecondMomentScratch[i][j];
                out.dLogL_dPrecisionXY[i][j] = -0.5 * crossSecondMomentScratch[j][i];
                out.dLogL_dPrecisionYX[i][j] = -0.5 * crossSecondMomentScratch[i][j];
                out.dLogL_dPrecisionYY[i][j] = -0.5 * nextSecondMomentScratch[i][j];
            }
        }
        out.dLogL_dLogNormalizer = -1.0;
    }

    private void multiplyMatrixMatrix(final double[][] left,
                                      final double[][] right,
                                      final double[][] out) {
        final int rows = left.length;
        final int inner = right.length;
        final int cols = right[0].length;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double sum = 0.0;
                for (int k = 0; k < inner; ++k) {
                    sum += left[i][k] * right[k][j];
                }
                out[i][j] = sum;
            }
        }
    }

    private void multiplyMatrixMatrixTransposedRight(final double[][] left,
                                                     final double[][] right,
                                                     final double[][] out) {
        final int rows = left.length;
        final int inner = left[0].length;
        final int cols = right.length;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double sum = 0.0;
                for (int k = 0; k < inner; ++k) {
                    sum += left[i][k] * right[j][k];
                }
                out[i][j] = sum;
            }
        }
    }

    private void subtractInPlace(final double[][] target, final double[][] delta) {
        for (int i = 0; i < target.length; ++i) {
            for (int j = 0; j < target[i].length; ++j) {
                target[i][j] -= delta[i][j];
            }
        }
    }

    private void addInPlace(final double[][] target, final double[][] delta) {
        for (int i = 0; i < target.length; ++i) {
            for (int j = 0; j < target[i].length; ++j) {
                target[i][j] += delta[i][j];
            }
        }
    }
}
