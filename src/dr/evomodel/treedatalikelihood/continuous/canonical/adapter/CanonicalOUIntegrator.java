/*
 * CanonicalOUIntegrator.java
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
 * Centralized backend for canonical OU tree likelihood state.
 *
 * <p>This is the canonical analogue of the continuous-diffusion integrator
 * ownership boundary: callers drive likelihood, gradients, tip refresh, and
 * MCMC state lifecycle through one object without depending on the concrete
 * message-passing implementation.</p>
 */
public interface CanonicalOUIntegrator {

    double calculateLogLikelihood();

    void computeSelectionGradient(double[] gradA);

    void computeDiffusionGradient(double[] gradQ);

    void computeStationaryMeanGradient(double[] gradMu);

    void computeGradientRootMean(double[] gradRootMean);

    void computeGradientBranchLengths(double[] gradT);

    void markTipObservationsDirty();

    void markModelDirty();

    void storeState();

    void restoreState();

    void acceptState();

}
