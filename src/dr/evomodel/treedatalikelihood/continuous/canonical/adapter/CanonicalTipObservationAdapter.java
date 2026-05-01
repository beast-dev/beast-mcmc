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

package dr.evomodel.treedatalikelihood.continuous.canonical.adapter;

import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalTipObservation;

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

    public static void fillTipObservations(final Tree tree,
                                           final ContinuousTraitPartialsProvider dataModel,
                                           final int dim,
                                           final CanonicalTipObservation[] observations) {
        if (dataModel.getTraitCount() != 1) {
            throw new UnsupportedOperationException(
                    "Canonical OU XML wiring currently supports only a single trait partition.");
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
                "Canonical OU XML wiring supports only exact observed-or-missing tips; "
                        + "encountered finite tip precision " + diagonalPrecision + " at index " + index + ".");
    }
}
