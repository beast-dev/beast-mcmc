/*
 * CanonicalTreeMessagePasser.java
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

package dr.evomodel.treedatalikelihood.continuous.framework;

import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianState;

/**
 * Canonical-form tree message passer contract.
 *
 * <p>The canonical node messages are represented directly by
 * {@link CanonicalGaussianState}. Exact tip observations are supplied through
 * {@link CanonicalTipObservation}, which keeps the canonical pathway free of the
 * variance-form conventions used by the older tree likelihood code.
 */
public interface CanonicalTreeMessagePasser {

    int getDimension();

    int getTipCount();

    void setTipObservation(int tipIndex, CanonicalTipObservation observation);

    double computePostOrderLogLikelihood(CanonicalBranchTransitionProvider transitionProvider,
                                         CanonicalRootPrior rootPrior);

    CanonicalGaussianState getPostOrderState(int nodeIndex);

    void computePreOrder(CanonicalBranchTransitionProvider transitionProvider,
                         CanonicalRootPrior rootPrior);

    CanonicalGaussianState getPreOrderState(int nodeIndex);

    /**
     * Legacy single-target gradient entry point. Production canonical OU callers
     * should prefer {@link #computeJointGradients} through the integrator so
     * shared branch adjoints are prepared once.
     */
    @Deprecated
    void computeGradientQ(CanonicalBranchTransitionProvider transitionProvider, double[] gradQ);

    void computeGradientBranchLengths(CanonicalBranchTransitionProvider transitionProvider, double[] gradT);

    /**
     * Legacy single-target gradient entry point. Production canonical OU callers
     * should prefer {@link #computeJointGradients} through the integrator so
     * shared branch adjoints are prepared once.
     */
    @Deprecated
    void computeGradientA(CanonicalBranchTransitionProvider transitionProvider, double[] gradA);

    /**
     * Legacy single-target gradient entry point. Production canonical OU callers
     * should prefer {@link #computeJointGradients} through the integrator so
     * shared branch adjoints are prepared once.
     */
    @Deprecated
    void computeGradientMu(CanonicalBranchTransitionProvider transitionProvider, double[] gradMu);

    default void computeJointGradients(final CanonicalBranchTransitionProvider transitionProvider,
                                       final double[] gradA,
                                       final double[] gradQ,
                                       final double[] gradMu) {
        computeGradientA(transitionProvider, gradA);
        computeGradientQ(transitionProvider, gradQ);
        computeGradientMu(transitionProvider, gradMu);
    }

    void storeState();

    void restoreState();

    void acceptState();
}
