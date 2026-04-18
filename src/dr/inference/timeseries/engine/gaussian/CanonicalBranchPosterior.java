package dr.inference.timeseries.engine.gaussian;

import dr.inference.timeseries.representation.CanonicalGaussianState;
import dr.inference.timeseries.representation.CanonicalGaussianUtils;

/**
 * Joint posterior statistics for a single branch pair {@code (x_t, x_{t+1})}.
 *
 * <p>The canonical form is the primary representation. Moment-form accessors are
 * derived from it on demand for branch-gradient calculations.
 *
 * <p>All internal matrix storage uses flat row-major arrays.
 */
final class CanonicalBranchPosterior {

    private final int stateDimension;
    private final CanonicalGaussianState pairState;
    private final double[] jointMean;
    private final double[][] jointCovariance;          // kept double[][] for CanonicalGaussianUtils bridge
    private final double[] currentMeanScratch;
    private final double[] nextMeanScratch;
    private final double[][] currentCovarianceScratch;   // double[][] for CanonicalGaussianUtils bridge
    private final double[][] crossCovarianceScratch;
    private final double[][] nextCovarianceScratch;
    private final double[] currentSecondMomentScratch;   // flat
    private final double[] crossSecondMomentScratch;     // flat
    private final double[] nextSecondMomentScratch;      // flat

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
        this.currentSecondMomentScratch = new double[stateDimension * stateDimension];
        this.crossSecondMomentScratch   = new double[stateDimension * stateDimension];
        this.nextSecondMomentScratch    = new double[stateDimension * stateDimension];
    }

    void fillFromMoments(final double[] currentMean,
                         final double[][] currentCovariance,
                         final double[] nextMean,
                         final double[][] nextCovariance,
                         final double[][] crossCovariance) {
        final int d = stateDimension;
        for (int i = 0; i < d; ++i) {
            jointMean[i]     = currentMean[i];
            jointMean[d + i] = nextMean[i];
        }
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                jointCovariance[i][j]         = currentCovariance[i][j];
                jointCovariance[i][d + j]     = crossCovariance[j][i];
                jointCovariance[d + i][j]     = crossCovariance[i][j];
                jointCovariance[d + i][d + j] = nextCovariance[i][j];
            }
        }
        CanonicalGaussianUtils.fillStateFromMoments(jointMean, jointCovariance, pairState);
    }

    CanonicalGaussianState getPairState() {
        return pairState;
    }

    void fillFromCanonicalPairState(final CanonicalGaussianState source) {
        final int dim2 = pairState.getDimension();
        System.arraycopy(source.information, 0, pairState.information, 0, dim2);
        System.arraycopy(source.precision,   0, pairState.precision,   0, dim2 * dim2);
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
            nextMeanOut[i]    = jointMean[d + i];
            for (int j = 0; j < d; ++j) {
                currentCovarianceOut[i][j] = jointCovariance[i][j];
                crossCovarianceOut[i][j]   = jointCovariance[d + i][j];
                nextCovarianceOut[i][j]    = jointCovariance[d + i][d + j];
            }
        }
    }

    void fillSufficientStatistics(final double[] currentMeanOut,
                                  final double[] nextMeanOut,
                                  final double[] currentSecondMomentOut,
                                  final double[] crossSecondMomentOut,
                                  final double[] nextSecondMomentOut) {
        CanonicalGaussianUtils.fillMomentsFromCanonical(pairState, jointMean, jointCovariance);
        final int d = stateDimension;
        for (int i = 0; i < d; ++i) {
            final double currentMeanI = jointMean[i];
            final double nextMeanI    = jointMean[d + i];
            currentMeanOut[i] = currentMeanI;
            nextMeanOut[i]    = nextMeanI;
            for (int j = 0; j < d; ++j) {
                final int ij = i * d + j;
                currentSecondMomentOut[ij] = jointCovariance[i][j] + currentMeanI * jointMean[j];
                crossSecondMomentOut[ij]   = jointCovariance[d + i][j] + nextMeanI * jointMean[j];
                nextSecondMomentOut[ij]    = jointCovariance[d + i][d + j] + nextMeanI * jointMean[d + j];
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
                residualCrossMomentOut[i][j]  = crossCovarianceScratch[i][j];
                residualSecondMomentOut[i][j] = nextCovarianceScratch[i][j];
            }
        }
        multiplyMatrixMatrix(transitionMatrix, currentCovarianceScratch, currentSecondMomentScratch);
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                residualCrossMomentOut[i][j] -= currentSecondMomentScratch[i * d + j];
                residualCrossMomentOut[i][j] += meanResidualOut[i] * currentMeanScratch[j];
            }
        }
        multiplyMatrixMatrixTransposedRight(transitionMatrix, crossCovarianceScratch, currentSecondMomentScratch);
        subtractInPlace(residualSecondMomentOut, currentSecondMomentScratch, d);
        multiplyMatrixMatrixTransposedRight(crossCovarianceScratch, transitionMatrix, currentSecondMomentScratch);
        subtractInPlace(residualSecondMomentOut, currentSecondMomentScratch, d);
        multiplyMatrixMatrix(transitionMatrix, currentCovarianceScratch, currentSecondMomentScratch);
        multiplyMatrixMatrixTransposedRight2(currentSecondMomentScratch, transitionMatrix, nextSecondMomentScratch, d);
        addInPlace(residualSecondMomentOut, nextSecondMomentScratch, d);
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
                final int ij = i * d + j;
                out.dLogL_dPrecisionXX[ij] = -0.5 * currentSecondMomentScratch[ij];
                out.dLogL_dPrecisionXY[ij] = -0.5 * crossSecondMomentScratch[j * d + i];
                out.dLogL_dPrecisionYX[ij] = -0.5 * crossSecondMomentScratch[ij];
                out.dLogL_dPrecisionYY[ij] = -0.5 * nextSecondMomentScratch[ij];
            }
        }
        out.dLogL_dLogNormalizer = -1.0;
    }

    // -----------------------------------------------------------------------
    // Private matrix helpers (mixed flat/ragged as needed)
    // -----------------------------------------------------------------------

    /** ragged × ragged → flat */
    private void multiplyMatrixMatrix(final double[][] left,
                                      final double[][] right,
                                      final double[] out) {
        final int d = stateDimension;
        for (int i = 0; i < d; ++i) {
            final int iOff = i * d;
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += left[i][k] * right[k][j];
                }
                out[iOff + j] = sum;
            }
        }
    }

    /** ragged × ragged^T → flat */
    private void multiplyMatrixMatrixTransposedRight(final double[][] left,
                                                     final double[][] right,
                                                     final double[] out) {
        final int d = stateDimension;
        for (int i = 0; i < d; ++i) {
            final int iOff = i * d;
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += left[i][k] * right[j][k];
                }
                out[iOff + j] = sum;
            }
        }
    }

    /** flat × ragged^T → flat */
    private void multiplyMatrixMatrixTransposedRight2(final double[] left,
                                                      final double[][] right,
                                                      final double[] out,
                                                      final int d) {
        for (int i = 0; i < d; ++i) {
            final int iOff = i * d;
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += left[iOff + k] * right[j][k];
                }
                out[iOff + j] = sum;
            }
        }
    }

    private void subtractInPlace(final double[][] target, final double[] delta, final int d) {
        for (int i = 0; i < d; ++i) {
            final int iOff = i * d;
            for (int j = 0; j < d; ++j) {
                target[i][j] -= delta[iOff + j];
            }
        }
    }

    private void addInPlace(final double[][] target, final double[] delta, final int d) {
        for (int i = 0; i < d; ++i) {
            final int iOff = i * d;
            for (int j = 0; j < d; ++j) {
                target[i][j] += delta[iOff + j];
            }
        }
    }
}
