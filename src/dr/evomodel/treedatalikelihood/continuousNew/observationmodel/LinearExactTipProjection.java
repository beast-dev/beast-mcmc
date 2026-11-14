package dr.evomodel.treedatalikelihood.continuousNew.observationmodel;

import dr.evomodel.continuous.canonical.math.MatrixOps;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalTransitionMomentProvider;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalLocalTransitionAdjoints;

import java.util.Arrays;

/**
 * Allocation-free projection and reverse-mode adjoint for exact linear tip
 * observations through a moment-form branch transition.
 */
public final class LinearExactTipProjection {

    private static final double LOG_TWO_PI = Math.log(2.0 * Math.PI);

    private final int dimension;
    private final double[] transitionMatrix;
    private final double[] transitionOffset;
    private final double[] covariance;
    private final double[] hOmega;
    private final double[] effectiveShifted;
    private final double[] sCompact;
    private final double[] sInverse;
    private final double[] forwardWork;
    private final double[] backwardWork;
    private final double[] w;
    private final double[] gMatrix;
    private final double[] sInvTimesG;
    private final double[] gInfoBar;
    private final double[] b;
    private final double[] gBar;
    private final double[] sBar;
    private final double[] uPrec;
    private final double[] sBarH;
    private final double[] hBar;
    private final double[] shiftedBar;

    public LinearExactTipProjection(final int dimension) {
        if (dimension < 1) {
            throw new IllegalArgumentException("dimension must be positive");
        }
        this.dimension = dimension;
        final int d2 = dimension * dimension;
        this.transitionMatrix = new double[d2];
        this.transitionOffset = new double[dimension];
        this.covariance = new double[d2];
        this.hOmega = new double[d2];
        this.effectiveShifted = new double[dimension];
        this.sCompact = new double[d2];
        this.sInverse = new double[d2];
        this.forwardWork = new double[dimension];
        this.backwardWork = new double[dimension];
        this.w = new double[dimension];
        this.gMatrix = new double[d2];
        this.sInvTimesG = new double[d2];
        this.gInfoBar = new double[dimension];
        this.b = new double[dimension];
        this.gBar = new double[d2];
        this.sBar = new double[d2];
        this.uPrec = new double[d2];
        this.sBarH = new double[d2];
        this.hBar = new double[d2];
        this.shiftedBar = new double[dimension];
    }

    public void projectObservedLinearChildToParent(
            final LinearExactTipObservationModel tipObservation,
            final int childNodeIndex,
            final CanonicalTransitionMomentProvider momentProvider,
            final double branchLength,
            final CanonicalGaussianState out) {
        final int observedCount = fillForwardTerms(
                tipObservation, childNodeIndex, momentProvider, branchLength);
        fillProjectedState(observedCount, out);
    }

    public void fillLocalAdjoints(
            final LinearExactTipObservationModel tipObservation,
            final int childNodeIndex,
            final CanonicalTransitionMomentProvider momentProvider,
            final double branchLength,
            final double[] informationBar,
            final double[] precisionBar,
            final double logNormalizerBar,
            final CanonicalLocalTransitionAdjoints out) {
        final int observedCount = fillForwardTerms(
                tipObservation, childNodeIndex, momentProvider, branchLength);
        fillAdjointsFromSeeds(
                tipObservation,
                observedCount,
                informationBar,
                precisionBar,
                logNormalizerBar,
                out);
    }

    public void fillLocalAdjointsForCurrentProjection(
            final LinearExactTipObservationModel tipObservation,
            final double[] informationBar,
            final double[] precisionBar,
            final double logNormalizerBar,
            final CanonicalLocalTransitionAdjoints out) {
        if (tipObservation.getLatentDimension() != dimension) {
            throw new IllegalArgumentException("Tip latent dimension mismatch");
        }
        fillAdjointsFromSeeds(
                tipObservation,
                tipObservation.getObservationDimension(),
                informationBar,
                precisionBar,
                logNormalizerBar,
                out);
    }

    private int fillForwardTerms(final LinearExactTipObservationModel tipObservation,
                                 final int childNodeIndex,
                                 final CanonicalTransitionMomentProvider momentProvider,
                                 final double branchLength) {
        if (tipObservation.getLatentDimension() != dimension) {
            throw new IllegalArgumentException("Tip latent dimension mismatch");
        }
        final int observedCount = tipObservation.getObservationDimension();
        momentProvider.fillTransitionMomentsFlat(
                childNodeIndex,
                branchLength,
                transitionMatrix,
                transitionOffset,
                covariance);

        for (int obs = 0; obs < observedCount; ++obs) {
            double offset = 0.0;
            final int obsOffset = obs * dimension;
            for (int j = 0; j < dimension; ++j) {
                offset += tipObservation.linkValueAt(obs, j) * transitionOffset[j];
            }
            effectiveShifted[obs] = tipObservation.shiftedValueAt(obs) - offset;

            for (int col = 0; col < dimension; ++col) {
                double g = 0.0;
                double hCov = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    final double h = tipObservation.linkValueAt(obs, k);
                    g += h * transitionMatrix[k * dimension + col];
                    hCov += h * covariance[k * dimension + col];
                }
                gMatrix[obsOffset + col] = g;
                hOmega[obsOffset + col] = hCov;
            }
        }

        for (int row = 0; row < observedCount; ++row) {
            final int rowOffset = row * observedCount;
            for (int col = 0; col < observedCount; ++col) {
                double sum = 0.0;
                final int hOmegaRow = row * dimension;
                for (int k = 0; k < dimension; ++k) {
                    sum += hOmega[hOmegaRow + k] * tipObservation.linkValueAt(col, k);
                }
                sCompact[rowOffset + col] = sum;
            }
        }
        symmetrizeCompact(sCompact, observedCount);

        final double logDet = MatrixOps.invertSPDCompact(
                sCompact,
                sInverse,
                observedCount,
                forwardWork,
                backwardWork);

        for (int row = 0; row < observedCount; ++row) {
            double sum = 0.0;
            final int rowOffset = row * observedCount;
            for (int k = 0; k < observedCount; ++k) {
                sum += sInverse[rowOffset + k] * effectiveShifted[k];
            }
            w[row] = sum;
        }

        for (int row = 0; row < observedCount; ++row) {
            final int rowOffset = row * dimension;
            final int sRowOffset = row * observedCount;
            for (int col = 0; col < dimension; ++col) {
                double sum = 0.0;
                for (int k = 0; k < observedCount; ++k) {
                    sum += sInverse[sRowOffset + k] * gMatrix[k * dimension + col];
                }
                sInvTimesG[rowOffset + col] = sum;
            }
        }

        backwardWork[0] = logDet;
        return observedCount;
    }

    private void fillProjectedState(final int observedCount,
                                    final CanonicalGaussianState out) {
        for (int col = 0; col < dimension; ++col) {
            double information = 0.0;
            for (int obs = 0; obs < observedCount; ++obs) {
                information += gMatrix[obs * dimension + col] * w[obs];
            }
            out.information[col] = information;
        }

        for (int i = 0; i < dimension; ++i) {
            final int iOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                double precision = 0.0;
                for (int obs = 0; obs < observedCount; ++obs) {
                    precision += gMatrix[obs * dimension + i]
                            * sInvTimesG[obs * dimension + j];
                }
                out.precision[iOffset + j] = precision;
            }
        }
        MatrixOps.symmetrize(out.precision, dimension);

        double quadratic = 0.0;
        for (int obs = 0; obs < observedCount; ++obs) {
            quadratic += effectiveShifted[obs] * w[obs];
        }
        out.logNormalizer = 0.5 * (observedCount * LOG_TWO_PI + backwardWork[0] + quadratic);
    }

    private void fillAdjointsFromSeeds(final LinearExactTipObservationModel tipObservation,
                                       final int observedCount,
                                       final double[] informationBar,
                                       final double[] precisionBar,
                                       final double logNormalizerBar,
                                       final CanonicalLocalTransitionAdjoints out) {
        Arrays.fill(out.dLogL_dF, 0.0);
        Arrays.fill(out.dLogL_df, 0.0);
        Arrays.fill(out.dLogL_dOmega, 0.0);
        out.clearObservationBasisAdjoint();

        for (int obs = 0; obs < observedCount; ++obs) {
            double sum = 0.0;
            final int obsOffset = obs * dimension;
            for (int j = 0; j < dimension; ++j) {
                sum += gMatrix[obsOffset + j] * informationBar[j];
            }
            gInfoBar[obs] = sum;
        }

        for (int obs = 0; obs < observedCount; ++obs) {
            double sum = 0.0;
            final int obsOffset = obs * observedCount;
            for (int j = 0; j < observedCount; ++j) {
                sum += sInverse[obsOffset + j] * gInfoBar[j];
            }
            b[obs] = sum;
        }

        for (int obs = 0; obs < observedCount; ++obs) {
            final int obsOffset = obs * dimension;
            for (int col = 0; col < dimension; ++col) {
                double sandwich = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sandwich += sInvTimesG[obsOffset + k] * precisionBar[k * dimension + col];
                }
                uPrec[obsOffset + col] = sandwich;
                gBar[obsOffset + col] = w[obs] * informationBar[col] + 2.0 * sandwich;
            }
        }

        for (int row = 0; row < observedCount; ++row) {
            final int rowOffset = row * observedCount;
            for (int col = 0; col < observedCount; ++col) {
                double value = -0.5 * (b[row] * w[col] + w[row] * b[col]);
                double sandwich = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sandwich += uPrec[row * dimension + k]
                            * sInvTimesG[col * dimension + k];
                }
                value -= sandwich;
                value += 0.5 * logNormalizerBar
                        * (sInverse[rowOffset + col] - w[row] * w[col]);
                sBar[rowOffset + col] = value;
            }
        }
        symmetrizeCompact(sBar, observedCount);

        for (int obs = 0; obs < observedCount; ++obs) {
            shiftedBar[obs] = b[obs] + logNormalizerBar * w[obs];
        }

        fillTransitionMatrixAdjoint(tipObservation, observedCount, out);
        fillCovarianceAdjoint(tipObservation, observedCount, out);
        fillObservationAndMeanAdjoints(tipObservation, observedCount, out);
        out.hasObservationBasisAdjoint = true;
    }

    private void fillTransitionMatrixAdjoint(final LinearExactTipObservationModel tipObservation,
                                             final int observedCount,
                                             final CanonicalLocalTransitionAdjoints out) {
        for (int row = 0; row < dimension; ++row) {
            final int rowOffset = row * dimension;
            for (int col = 0; col < dimension; ++col) {
                double sum = 0.0;
                for (int obs = 0; obs < observedCount; ++obs) {
                    sum += tipObservation.linkValueAt(obs, row)
                            * gBar[obs * dimension + col];
                }
                out.dLogL_dF[rowOffset + col] = sum;
            }
        }
    }

    private void fillCovarianceAdjoint(final LinearExactTipObservationModel tipObservation,
                                       final int observedCount,
                                       final CanonicalLocalTransitionAdjoints out) {
        for (int obs = 0; obs < observedCount; ++obs) {
            final int obsOffset = obs * dimension;
            final int sOffset = obs * observedCount;
            for (int col = 0; col < dimension; ++col) {
                double sum = 0.0;
                for (int k = 0; k < observedCount; ++k) {
                    sum += sBar[sOffset + k] * tipObservation.linkValueAt(k, col);
                }
                sBarH[obsOffset + col] = sum;
            }
        }

        for (int row = 0; row < dimension; ++row) {
            final int rowOffset = row * dimension;
            for (int col = 0; col < dimension; ++col) {
                double sum = 0.0;
                for (int obs = 0; obs < observedCount; ++obs) {
                    sum += tipObservation.linkValueAt(obs, row)
                            * sBarH[obs * dimension + col];
                }
                out.dLogL_dOmega[rowOffset + col] = sum;
            }
        }
        MatrixOps.symmetrize(out.dLogL_dOmega, dimension);
    }

    private void fillObservationAndMeanAdjoints(final LinearExactTipObservationModel tipObservation,
                                                final int observedCount,
                                                final CanonicalLocalTransitionAdjoints out) {
        for (int obs = 0; obs < observedCount; ++obs) {
            final int obsOffset = obs * dimension;
            for (int col = 0; col < dimension; ++col) {
                double fromG = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    fromG += gBar[obsOffset + k] * transitionMatrix[col * dimension + k];
                }

                double fromS = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    fromS += sBarH[obsOffset + k] * covariance[k * dimension + col];
                }

                hBar[obsOffset + col] =
                        fromG + 2.0 * fromS - shiftedBar[obs] * transitionOffset[col];
            }
        }

        for (int obs = 0; obs < observedCount; ++obs) {
            final int naturalIndex = tipObservation.observedIndexAt(obs);
            final int rOffset = naturalIndex * dimension;
            for (int col = 0; col < dimension; ++col) {
                out.dLogL_dObservationBasis[rOffset + col] += hBar[obs * dimension + col];
            }
            out.dLogL_dStationaryMeanDirect[naturalIndex] -= shiftedBar[obs];
        }

        for (int col = 0; col < dimension; ++col) {
            double offsetAdjoint = 0.0;
            for (int obs = 0; obs < observedCount; ++obs) {
                offsetAdjoint -= tipObservation.linkValueAt(obs, col) * shiftedBar[obs];
            }
            out.dLogL_df[col] = offsetAdjoint;
        }
    }

    private static void symmetrizeCompact(final double[] matrix, final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = i + 1; j < dimension; ++j) {
                final double value = 0.5 * (matrix[i * dimension + j] + matrix[j * dimension + i]);
                matrix[i * dimension + j] = value;
                matrix[j * dimension + i] = value;
            }
        }
    }
}
