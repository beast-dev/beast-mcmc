/*
 * CanonicalTransitionCacheEpoch.java
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

/**
 * Version counters used to validate branch transition-cache entries.
 */
final class CanonicalTransitionCacheEpoch {

    private long selectionVersion;
    private long diffusionVersion;
    private long stationaryMeanVersion;
    private long branchLengthVersion;

    long selectionVersion() {
        return selectionVersion;
    }

    long diffusionVersion() {
        return diffusionVersion;
    }

    long stationaryMeanVersion() {
        return stationaryMeanVersion;
    }

    long branchLengthVersion() {
        return branchLengthVersion;
    }

    void advance(final CanonicalTransitionCacheInvalidationReason reason) {
        if (reason == CanonicalTransitionCacheInvalidationReason.DIFFUSION_CHANGED) {
            diffusionVersion++;
        } else if (reason == CanonicalTransitionCacheInvalidationReason.STATIONARY_MEAN_CHANGED) {
            stationaryMeanVersion++;
        } else if (reason == CanonicalTransitionCacheInvalidationReason.SELECTION_CHANGED) {
            selectionVersion++;
        } else if (reason == CanonicalTransitionCacheInvalidationReason.BRANCH_LENGTH_CHANGED) {
            branchLengthVersion++;
        } else {
            selectionVersion++;
            diffusionVersion++;
            stationaryMeanVersion++;
            branchLengthVersion++;
        }
    }
}
