/*
 * DiscontinuousPotentialProvider.java
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

package dr.inference.hmc;

import dr.inference.model.Parameter;
import dr.inference.operators.RewardMixturePerformanceStats;

/**
 * Provides local log-density information needed by discontinuous HMC updates.
 *
 * The initial implementation focuses on coordinatewise discontinuous updates, so
 * the key operation is evaluation of the log density after moving a single
 * coordinate while holding the others fixed.
 *
 * @author Filippo Monti (powered by OpenAI)
 */
public interface DiscontinuousPotentialProvider {

    Parameter getParameter();

    int getDimension();

    /**
     * Called exactly once per HMC operator call, before any coordinate
     * integration, by callers that own the top of that call (e.g.
     * {@code MixedDiscontinuousHamiltonianMonteCarloOperator.doOperation}).
     * Providers whose potential landscape depends on auxiliary state that
     * must stay static for the duration of one trajectory (see
     * {@code PerBranchRewardMixtureCategoryDecoder}) should snapshot that
     * state here rather than inside per-crossing methods like
     * {@link #getPotentialDifferenceAcrossBoundary}, which are called many
     * times per operator call and must not repeat O(problem size) work.
     */
    default void refresh() {
    }

    /**
     * Called after smooth position integration has written the final proposal
     * coordinates, immediately before final proposal scoring. Most
     * discontinuous providers have no final-position cache to synchronize.
     */
    default void refreshAfterPositionUpdate() {
    }

    /**
     * Called after an integrated position update has changed model state while
     * the current operation is still in progress. Providers can drop derived
     * likelihood/message caches here without changing any frozen discontinuity
     * layout captured by {@link #refresh()}.
     */
    default void clearOperationCache() {
    }

    default void clearOperationCache(final RewardMixturePerformanceStats.OperationCacheClearReason reason) {
        clearOperationCache();
    }

    double getLogDensity();

    double getLogDensityAfterSingleCoordinateMove(int index, double proposedValue);

    default boolean isDiscontinuous(int index) {
        return true;
    }

    default double getNextDiscontinuity(int index, double currentValue, double direction) {
        throw new UnsupportedOperationException(
                getClass().getName() + " does not expose discontinuity locations");
    }

    default double getPotentialDifferenceAcrossBoundary(int index, double boundary, double direction) {
        if (direction == 0.0) {
            return 0.0;
        }

        final double before = Math.nextAfter(boundary, boundary - direction);
        final double after = Math.nextAfter(boundary, boundary + direction);
        final double currentLogDensity = getLogDensityAfterSingleCoordinateMove(index, before);
        final double proposedLogDensity = getLogDensityAfterSingleCoordinateMove(index, after);

        if (Double.isInfinite(proposedLogDensity) && proposedLogDensity < 0.0) {
            return Double.POSITIVE_INFINITY;
        }
        if (Double.isInfinite(currentLogDensity) && currentLogDensity < 0.0) {
            return Double.NEGATIVE_INFINITY;
        }

        return currentLogDensity - proposedLogDensity;
    }

    default double getPotentialDifference(int index, double currentValue, double proposedValue) {
        final double currentLogDensity = getLogDensity();
        final double proposedLogDensity = getLogDensityAfterSingleCoordinateMove(index, proposedValue);

        if (Double.isInfinite(proposedLogDensity) && proposedLogDensity < 0.0) {
            return Double.POSITIVE_INFINITY;
        }

        return currentLogDensity - proposedLogDensity;
    }
}
