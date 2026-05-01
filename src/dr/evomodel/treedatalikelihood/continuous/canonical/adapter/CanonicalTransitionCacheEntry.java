/*
 * CanonicalTransitionCacheEntry.java
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

import dr.evomodel.continuous.ou.canonical.CanonicalPreparedBranchHandle;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalPreparedBranchSnapshot;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;

/**
 * Mutable per-branch entry owned by {@link CanonicalTransitionCache}.
 */
final class CanonicalTransitionCacheEntry {

    CanonicalGaussianTransition transition;
    CanonicalPreparedBranchSnapshot snapshot;
    CanonicalPreparedBranchHandle preparedBranchHandle;
    double effectiveBranchLength;
    boolean preparedBasisValid;

    private boolean valid;
    private long transitionSelectionVersion;
    private long transitionDiffusionVersion;
    private long transitionStationaryMeanVersion;
    private long transitionBranchLengthVersion;
    private long preparedSelectionVersion;
    private long preparedStationaryMeanVersion;
    private long preparedBranchLengthVersion;

    boolean isTransitionCurrent(final double effectiveBranchLength,
                                final CanonicalTransitionCacheEpoch epoch) {
        return valid
                && Double.doubleToLongBits(this.effectiveBranchLength)
                == Double.doubleToLongBits(effectiveBranchLength)
                && transitionSelectionVersion == epoch.selectionVersion()
                && transitionDiffusionVersion == epoch.diffusionVersion()
                && transitionStationaryMeanVersion == epoch.stationaryMeanVersion()
                && transitionBranchLengthVersion == epoch.branchLengthVersion();
    }

    boolean hasBranchLengthChanged(final double effectiveBranchLength) {
        return valid
                && Double.doubleToLongBits(this.effectiveBranchLength)
                != Double.doubleToLongBits(effectiveBranchLength);
    }

    void markTransitionCurrent(final double effectiveBranchLength,
                               final CanonicalTransitionCacheEpoch epoch) {
        this.effectiveBranchLength = effectiveBranchLength;
        this.transitionSelectionVersion = epoch.selectionVersion();
        this.transitionDiffusionVersion = epoch.diffusionVersion();
        this.transitionStationaryMeanVersion = epoch.stationaryMeanVersion();
        this.transitionBranchLengthVersion = epoch.branchLengthVersion();
        this.valid = true;
    }

    boolean isPreparedBasisCurrent(final CanonicalTransitionCacheEpoch epoch) {
        return preparedBasisValid
                && preparedSelectionVersion == epoch.selectionVersion()
                && preparedStationaryMeanVersion == epoch.stationaryMeanVersion()
                && preparedBranchLengthVersion == epoch.branchLengthVersion();
    }

    void markPreparedBasisCurrent(final CanonicalTransitionCacheEpoch epoch) {
        this.preparedSelectionVersion = epoch.selectionVersion();
        this.preparedStationaryMeanVersion = epoch.stationaryMeanVersion();
        this.preparedBranchLengthVersion = epoch.branchLengthVersion();
        this.preparedBasisValid = true;
    }

    boolean invalidate(final CanonicalTransitionCacheInvalidationReason reason) {
        valid = false;
        if (preparedBranchHandle == null) {
            preparedBasisValid = false;
            return false;
        }
        if (reason == CanonicalTransitionCacheInvalidationReason.DIFFUSION_CHANGED) {
            preparedBranchHandle.invalidateCovariance();
            return true;
        }
        preparedBasisValid = false;
        preparedBranchHandle.invalidateCovariance();
        return true;
    }
}
