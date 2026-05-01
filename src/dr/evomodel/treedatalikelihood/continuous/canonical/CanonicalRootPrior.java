/*
 * CanonicalRootPrior.java
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

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;

/**
 * Canonical-form analogue of {@link RootPrior}.
 *
 * <p>The canonical OU path uses this interface to seed the root message for the
 * downward pass and to integrate the canonical upward root message into the final
 * log-likelihood without forcing the core passer API through a moment-form seam.
 */
public interface CanonicalRootPrior {

    int getDimension();

    void fillRootPriorState(double[][] traitCovariance, CanonicalGaussianState out);

    double computeLogMarginalLikelihood(CanonicalGaussianState rootMessage, double[][] traitCovariance);

    default boolean isFixedRoot() {
        return false;
    }

    default void fillFixedRootValue(double[] out) {
        throw new UnsupportedOperationException("This canonical root prior does not expose a fixed root value.");
    }

    default void accumulateRootMeanGradient(CanonicalGaussianState rootMessage,
                                           double[][] traitCovariance,
                                           double[] gradOut) {
        // No root-mean contribution by default.
    }

    default double getDiffusionScale() {
        return 0.0;
    }
}
