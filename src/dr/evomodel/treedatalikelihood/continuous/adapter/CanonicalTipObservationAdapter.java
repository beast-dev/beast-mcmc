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

package dr.evomodel.treedatalikelihood.continuous.adapter;

import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalTipObservation;

/**
 * Converts BEAST tip partials into exact observed-or-missing canonical tips.
 *
 * <p>This adapter intentionally only supports the exact-observation boundary used
 * by the canonical OU tree passer: each coordinate must be either observed
 * exactly (precision {@code +Inf}) or missing (precision {@code 0}).
 */
public final class CanonicalTipObservationAdapter {

    private CanonicalTipObservationAdapter() { }

    public static CanonicalTipObservation[] extractTipObservations(final Tree tree,
                                                                   final ContinuousTraitPartialsProvider dataModel,
                                                                   final int dim) {
        if (dataModel.getTraitCount() != 1) {
            throw new UnsupportedOperationException(
                    "Canonical OU XML wiring currently supports only a single trait partition.");
        }

        final int tipCount = tree.getExternalNodeCount();
        final CanonicalTipObservation[] observations = new CanonicalTipObservation[tipCount];
        final PrecisionType precisionType = dataModel.getPrecisionType();
        for (int tipIdx = 0; tipIdx < tipCount; tipIdx++) {
            final double[] rawPartial = dataModel.getTipPartial(tipIdx, false);
            observations[tipIdx] = extractObservation(rawPartial, precisionType, dim);
        }
        return observations;
    }

    static CanonicalTipObservation extractObservation(final double[] raw,
                                                      final PrecisionType precisionType,
                                                      final int dim) {
        final CanonicalTipObservation observation = new CanonicalTipObservation(dim);
        final boolean[] observedMask = new boolean[dim];

        for (int i = 0; i < dim; i++) {
            observation.values[i] = raw[i];
            observedMask[i] = isExactlyObserved(raw, precisionType, dim, i);
        }
        observation.setPartiallyObserved(observation.values, observedMask);
        return observation;
    }

    private static boolean isExactlyObserved(final double[] raw,
                                             final PrecisionType precisionType,
                                             final int dim,
                                             final int index) {
        final double diagonalPrecision;
        switch (precisionType) {
            case FULL:
                diagonalPrecision = raw[precisionType.getPrecisionOffset(dim) + index * dim + index];
                break;
            case MIXED:
                diagonalPrecision = raw[precisionType.getPrecisionOffset(dim) + index];
                break;
            case SCALAR:
                diagonalPrecision = raw[precisionType.getPrecisionOffset(dim)];
                break;
            default:
                throw new UnsupportedOperationException(
                        "Canonical OU XML wiring does not support precision type " + precisionType.getTag());
        }

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
