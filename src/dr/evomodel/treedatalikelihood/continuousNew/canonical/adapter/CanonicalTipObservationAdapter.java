/*
 * CanonicalTipObservationAdapter.java
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

package dr.evomodel.treedatalikelihood.continuousNew.canonical.adapter;

import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodel.treedatalikelihood.continuousNew.observationmodel.CanonicalTipObservation;
import dr.evomodel.treedatalikelihood.continuousNew.observationmodel.CanonicalTipObservationModel;
import dr.evomodel.treedatalikelihood.continuousNew.observationmodel.GaussianCanonicalTipObservation;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.Arrays;

/**
 * Converts BEAST tip partials into exact observed-or-missing canonical tips.
 *
 * <p>This adapter intentionally only supports the exact-observation boundary used
 * by the canonical OU tree passer: each coordinate must be either observed
 * exactly (precision {@code +Inf}) or missing (precision {@code 0}).
 */
public final class CanonicalTipObservationAdapter {

    private CanonicalTipObservationAdapter() { }

    /**
     * True if at least one tip carries observation error, i.e. a finite positive diagonal precision (as supplied by a
     * repeated-measures data model with a sampling precision).  Only FULL precision partials of a single trait
     * partition can carry observation error; any other data model reports false.
     */
    public static boolean hasObservationError(final Tree tree,
                                              final ContinuousTraitPartialsProvider dataModel,
                                              final int dim) {
        if (tree == null || dataModel == null || dataModel.getTraitCount() != 1) {
            return false;
        }
        final PrecisionType precisionType = dataModel.getPrecisionType();
        if (precisionType != PrecisionType.FULL) {
            return false;
        }
        final int precisionOffset = precisionType.getPrecisionOffset(dim);
        final int tipCount = tree.getExternalNodeCount();
        for (int tipIdx = 0; tipIdx < tipCount; tipIdx++) {
            final double[] raw = dataModel.getTipPartial(tipIdx, false);
            for (int i = 0, d = precisionOffset; i < dim; i++, d += dim + 1) {
                final double precision = raw[d];
                if (!Double.isInfinite(precision) && precision > 0.0) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Like {@link #fillTipObservations(Tree, ContinuousTraitPartialsProvider, int, CanonicalTipObservation[])},
     * but additionally accepts tips observed with finite Gaussian error (finite diagonal precision, as supplied by a
     * repeated-measures data model with a sampling precision).  Such a tip is written to {@code noisyModels[tip]} as a
     * {@link GaussianCanonicalTipObservation} (identity link restricted to the observed coordinates, covariance equal
     * to the inverse of the observed precision block) and recorded as missing in {@code observations[tip]}; every
     * other entry of {@code noisyModels} is set to {@code null}.  A tip may not mix exact and noisy coordinates.
     * Only FULL precision partials are supported.
     *
     * @return true if at least one tip carries observation error
     */
    public static boolean fillTipObservations(final Tree tree,
                                              final ContinuousTraitPartialsProvider dataModel,
                                              final int dim,
                                              final CanonicalTipObservation[] observations,
                                              final CanonicalTipObservationModel[] noisyModels) {
        if (noisyModels == null) {
            fillTipObservations(tree, dataModel, dim, observations);
            return false;
        }
        if (dataModel.getTraitCount() != 1) {
            throw new UnsupportedOperationException(
                    "Canonical OU XML wiring requires a single trait partition.");
        }
        final int tipCount = tree.getExternalNodeCount();
        if (observations.length != tipCount || noisyModels.length != tipCount) {
            throw new IllegalArgumentException("Tip observation count mismatch");
        }
        final PrecisionType precisionType = dataModel.getPrecisionType();
        if (precisionType != PrecisionType.FULL) {
            fillTipObservations(tree, dataModel, dim, observations);
            Arrays.fill(noisyModels, null);
            return false;
        }
        final int precisionOffset = precisionType.getPrecisionOffset(dim);
        boolean anyNoisy = false;
        for (int tipIdx = 0; tipIdx < tipCount; tipIdx++) {
            final CanonicalTipObservation observation = validateObservation(observations, tipIdx, dim);
            final double[] raw = dataModel.getTipPartial(tipIdx, false);
            noisyModels[tipIdx] = null;
            int exact = 0;
            int noisy = 0;
            for (int i = 0, d = precisionOffset; i < dim; i++, d += dim + 1) {
                final double precision = raw[d];
                if (Double.isInfinite(precision)) {
                    exact++;
                } else if (precision > 0.0) {
                    noisy++;
                }
            }
            if (noisy == 0) {
                fillObservationFull(raw, dim, precisionOffset, observation);
                continue;
            }
            if (exact > 0) {
                throw new UnsupportedOperationException(
                        "Canonical OU tip " + tipIdx + " mixes exact and noisy coordinates.");
            }
            final int[] observedIndex = new int[noisy];
            int k = 0;
            for (int i = 0, d = precisionOffset; i < dim; i++, d += dim + 1) {
                if (raw[d] > 0.0) {
                    observedIndex[k++] = i;
                }
            }
            final DenseMatrix64F precisionBlock = new DenseMatrix64F(noisy, noisy);
            final double[] values = new double[noisy];
            final double[] link = new double[noisy * dim];
            for (int a = 0; a < noisy; a++) {
                values[a] = raw[observedIndex[a]];
                link[a * dim + observedIndex[a]] = 1.0;
                for (int b = 0; b < noisy; b++) {
                    precisionBlock.set(a, b, raw[precisionOffset + observedIndex[a] * dim + observedIndex[b]]);
                }
            }
            final DenseMatrix64F covariance = new DenseMatrix64F(noisy, noisy);
            if (!CommonOps.invert(precisionBlock, covariance)) {
                throw new IllegalStateException("Singular tip precision block at tip " + tipIdx);
            }
            noisyModels[tipIdx] = new GaussianCanonicalTipObservation(
                    dim, noisy, values, link, new double[noisy], covariance.data);
            observation.setMissing();
            anyNoisy = true;
        }
        return anyNoisy;
    }

    public static void fillTipObservations(final Tree tree,
                                           final ContinuousTraitPartialsProvider dataModel,
                                           final int dim,
                                           final CanonicalTipObservation[] observations) {
        if (dataModel.getTraitCount() != 1) {
            throw new UnsupportedOperationException(
                    "Canonical OU XML wiring requires a single trait partition.");
        }

        final int tipCount = tree.getExternalNodeCount();
        if (observations.length != tipCount) {
            throw new IllegalArgumentException(
                    "Tip observation count mismatch: " + observations.length + " vs " + tipCount);
        }
        final PrecisionType precisionType = dataModel.getPrecisionType();
        final int precisionOffset = precisionType.getPrecisionOffset(dim);

        switch (precisionType) {
            case FULL:
                for (int tipIdx = 0; tipIdx < tipCount; tipIdx++) {
                    final CanonicalTipObservation observation = validateObservation(observations, tipIdx, dim);
                    fillObservationFull(dataModel.getTipPartial(tipIdx, false), dim, precisionOffset, observation);
                }
                break;
            case MIXED:
                for (int tipIdx = 0; tipIdx < tipCount; tipIdx++) {
                    final CanonicalTipObservation observation = validateObservation(observations, tipIdx, dim);
                    fillObservationMixed(dataModel.getTipPartial(tipIdx, false), dim, precisionOffset, observation);
                }
                break;
            case SCALAR:
                for (int tipIdx = 0; tipIdx < tipCount; tipIdx++) {
                    final CanonicalTipObservation observation = validateObservation(observations, tipIdx, dim);
                    fillObservationScalar(dataModel.getTipPartial(tipIdx, false), dim, precisionOffset, observation);
                }
                break;
            default:
                throw new UnsupportedOperationException(
                        "Canonical OU XML wiring does not support precision type " + precisionType.getTag());
        }
    }

    static void fillObservation(final double[] raw,
                                final PrecisionType precisionType,
                                final int dim,
                                final CanonicalTipObservation observation) {
        final int precisionOffset = precisionType.getPrecisionOffset(dim);
        switch (precisionType) {
            case FULL:
                fillObservationFull(raw, dim, precisionOffset, observation);
                return;
            case MIXED:
                fillObservationMixed(raw, dim, precisionOffset, observation);
                return;
            case SCALAR:
                fillObservationScalar(raw, dim, precisionOffset, observation);
                return;
            default:
                throw new UnsupportedOperationException(
                        "Canonical OU XML wiring does not support precision type " + precisionType.getTag());
        }
    }

    private static CanonicalTipObservation validateObservation(final CanonicalTipObservation[] observations,
                                                               final int tipIdx,
                                                               final int dim) {
        final CanonicalTipObservation observation = observations[tipIdx];
        if (observation == null) {
            throw new IllegalArgumentException("observations[" + tipIdx + "] must not be null");
        }
        if (observation.dim != dim) {
            throw new IllegalArgumentException(
                    "Observation dimension mismatch at tip " + tipIdx + ": " + observation.dim + " vs " + dim);
        }
        return observation;
    }

    private static void fillObservationFull(final double[] raw,
                                            final int dim,
                                            final int precisionOffset,
                                            final CanonicalTipObservation observation) {
        System.arraycopy(raw, 0, observation.values, 0, dim);
        int observedCount = 0;
        int diagonalIndex = precisionOffset;
        for (int i = 0; i < dim; i++, diagonalIndex += dim + 1) {
            final boolean isObserved = decodeExactObservation(raw[diagonalIndex], i);
            observation.observed[i] = isObserved;
            if (isObserved) {
                observedCount++;
            }
        }
        observation.observedCount = observedCount;
    }

    private static void fillObservationMixed(final double[] raw,
                                             final int dim,
                                             final int precisionOffset,
                                             final CanonicalTipObservation observation) {
        System.arraycopy(raw, 0, observation.values, 0, dim);
        int observedCount = 0;
        for (int i = 0; i < dim; i++) {
            final boolean isObserved = decodeExactObservation(raw[precisionOffset + i], i);
            observation.observed[i] = isObserved;
            if (isObserved) {
                observedCount++;
            }
        }
        observation.observedCount = observedCount;
    }

    private static void fillObservationScalar(final double[] raw,
                                              final int dim,
                                              final int precisionOffset,
                                              final CanonicalTipObservation observation) {
        System.arraycopy(raw, 0, observation.values, 0, dim);
        final boolean isObserved = decodeExactObservation(raw[precisionOffset], 0);
        Arrays.fill(observation.observed, isObserved);
        observation.observedCount = isObserved ? dim : 0;
    }

    private static boolean decodeExactObservation(final double diagonalPrecision,
                                                  final int index) {
        if (Double.isInfinite(diagonalPrecision)) {
            return true;
        }
        if (diagonalPrecision == 0.0) {
            return false;
        }
        throw new UnsupportedOperationException(
                "Canonical OU XML wiring supports only exact observed-or-missing tips here; "
                        + "tip observation error needs canonicalOUCoordinates=\"natural\" and a repeatedMeasuresModel. "
                        + "Encountered finite tip precision " + diagonalPrecision + " at index " + index + ".");
    }
}
