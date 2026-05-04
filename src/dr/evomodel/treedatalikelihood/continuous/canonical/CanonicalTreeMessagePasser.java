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

package dr.evomodel.treedatalikelihood.continuous.canonical;

import dr.evomodel.treedatalikelihood.continuous.canonical.gradient.BranchGradientInputs;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.CanonicalTipObservation;

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

    void computeGradientBranchLengths(CanonicalBranchTransitionProvider transitionProvider, double[] gradT);

    void computeGradientBranchLengths(CanonicalBranchTransitionProvider transitionProvider,
                                      BranchGradientInputs inputs,
                                      double[] gradT);

    void computeJointGradients(CanonicalBranchTransitionProvider transitionProvider,
                               double[] gradA,
                               double[] gradQ,
                               double[] gradMu);

    BranchGradientInputs createBranchGradientInputs();

    void prepareBranchGradientInputs(CanonicalBranchTransitionProvider transitionProvider,
                                     String phase,
                                     BranchGradientInputs out);

    void computeJointGradients(CanonicalBranchTransitionProvider transitionProvider,
                               BranchGradientInputs inputs,
                               double[] gradA,
                               double[] gradQ,
                               double[] gradMu);

    void storeState();

    void restoreState();

    void acceptState();
}
