/*
 * CanonicalTipProjector.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.evomodel.treedatalikelihood.continuous.canonical;

import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalTipObservation;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalTransitionMomentProvider;
import dr.evomodel.treedatalikelihood.continuous.canonical.MatrixUtils;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianMessageOps;

/**
 * Projects exact tip observations through an OU branch into a parent-side
 * canonical message.
 */
final class CanonicalTipProjector {

    private static final double LOG_TWO_PI = Math.log(2.0 * Math.PI);

    private final int dimension;
    private final int[] observedIndexScratch;
    private final double[] shiftedObservationScratch;
    private final double[] precisionTimesShiftedScratch;
    private final double[] varianceFlatScratch;
    private final double[] precisionFlatScratch;
    private final double[] projectedPrecisionScratch;
    private final double[] choleskyScratch;
    private final double[] transitionMatrixFlat;
    private final double[] covarianceFlat;
    private final double[] transitionOffset;

    CanonicalTipProjector(final int dimension) {
        this.dimension = dimension;
        this.observedIndexScratch = new int[dimension];
        this.shiftedObservationScratch = new double[dimension];
        this.precisionTimesShiftedScratch = new double[dimension];
        this.varianceFlatScratch = new double[dimension * dimension];
        this.precisionFlatScratch = new double[dimension * dimension];
        this.projectedPrecisionScratch = new double[dimension * dimension];
        this.choleskyScratch = new double[dimension * dimension];
        this.transitionMatrixFlat = new double[dimension * dimension];
        this.covarianceFlat = new double[dimension * dimension];
        this.transitionOffset = new double[dimension];
    }

    void projectObservedChildToParent(final CanonicalTipObservation tipObservation,
                                      final CanonicalTransitionMomentProvider transitionMomentProvider,
                                      final double branchLength,
                                      final CanonicalGaussianState out) {
        final int observedCount = collectObservedIndices(tipObservation);
        if (observedCount == 0) {
            CanonicalGaussianMessageOps.clearState(out);
            return;
        }

        transitionMomentProvider.fillTransitionMatrixFlat(branchLength, transitionMatrixFlat);
        transitionMomentProvider.fillTransitionOffset(branchLength, transitionOffset);
        transitionMomentProvider.fillTransitionCovarianceFlat(branchLength, covarianceFlat);

        for (int observed = 0; observed < observedCount; ++observed) {
            final int observedTrait = observedIndexScratch[observed];
            shiftedObservationScratch[observed] =
                    tipObservation.values[observedTrait] - transitionOffset[observedTrait];
            final int rowOffset = observed * observedCount;
            for (int otherObserved = 0; otherObserved < observedCount; ++otherObserved) {
                varianceFlatScratch[rowOffset + otherObserved] =
                        covarianceFlat[observedTrait * dimension + observedIndexScratch[otherObserved]];
            }
        }

        final double logDetVariance = MatrixUtils.invertSymmetricPositiveDefiniteCompact(
                varianceFlatScratch,
                precisionFlatScratch,
                observedCount,
                transitionOffset,
                choleskyScratch);

        for (int row = 0; row < observedCount; ++row) {
            double sum = 0.0;
            final int rowOffset = row * observedCount;
            for (int k = 0; k < observedCount; ++k) {
                sum += precisionFlatScratch[rowOffset + k] * shiftedObservationScratch[k];
            }
            precisionTimesShiftedScratch[row] = sum;
        }

        for (int row = 0; row < observedCount; ++row) {
            final int precisionRowOffset = row * observedCount;
            final int projectedRowOffset = row * dimension;
            for (int col = 0; col < dimension; ++col) {
                double sum = 0.0;
                for (int k = 0; k < observedCount; ++k) {
                    sum += precisionFlatScratch[precisionRowOffset + k]
                            * transitionMatrixFlat[observedIndexScratch[k] * dimension + col];
                }
                projectedPrecisionScratch[projectedRowOffset + col] = sum;
            }
        }

        for (int i = 0; i < dimension; ++i) {
            double information = 0.0;
            for (int observed = 0; observed < observedCount; ++observed) {
                information += transitionMatrixFlat[observedIndexScratch[observed] * dimension + i]
                        * precisionTimesShiftedScratch[observed];
            }
            out.information[i] = information;

            for (int j = 0; j < dimension; ++j) {
                double precision = 0.0;
                for (int observed = 0; observed < observedCount; ++observed) {
                    precision += transitionMatrixFlat[observedIndexScratch[observed] * dimension + i]
                            * projectedPrecisionScratch[observed * dimension + j];
                }
                out.precision[i * dimension + j] = precision;
            }
        }
        symmetrizeFlatSquare(out.precision);

        double quadratic = 0.0;
        for (int observed = 0; observed < observedCount; ++observed) {
            quadratic += shiftedObservationScratch[observed] * precisionTimesShiftedScratch[observed];
        }
        out.logNormalizer = 0.5 * (observedCount * LOG_TWO_PI + logDetVariance + quadratic);
    }

    private int collectObservedIndices(final CanonicalTipObservation tipObservation) {
        int observedCount = 0;
        for (int i = 0; i < dimension; i++) {
            if (tipObservation.observed[i]) {
                observedIndexScratch[observedCount++] = i;
            }
        }
        if (observedCount != tipObservation.observedCount) {
            throw new UnsupportedOperationException(
                    "Canonical tip observation partition is inconsistent with observedCount.");
        }
        return observedCount;
    }

    private void symmetrizeFlatSquare(final double[] matrix) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = i + 1; j < dimension; ++j) {
                final double avg = 0.5 * (matrix[i * dimension + j] + matrix[j * dimension + i]);
                matrix[i * dimension + j] = avg;
                matrix[j * dimension + i] = avg;
            }
        }
    }
}
