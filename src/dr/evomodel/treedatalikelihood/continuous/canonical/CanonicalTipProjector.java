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
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.PartialIdentityTipProjection;

/**
 * Projects exact tip observations through an OU branch into a parent-side
 * canonical message.
 */
final class CanonicalTipProjector {

    private final PartialIdentityTipProjection partialIdentityProjection;

    CanonicalTipProjector(final int dimension) {
        this.partialIdentityProjection = new PartialIdentityTipProjection(dimension);
    }

    void projectObservedChildToParent(final CanonicalTipObservation tipObservation,
                                      final CanonicalTransitionMomentProvider transitionMomentProvider,
                                      final double branchLength,
                                      final CanonicalGaussianState out) {
        partialIdentityProjection.projectObservedChildToParent(
                tipObservation,
                transitionMomentProvider,
                branchLength,
                out);
    }

}
