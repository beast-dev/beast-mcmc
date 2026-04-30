/*
 * CanonicalOUGradientAdapter.java
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

import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;

/**
 * HMC-facing adapter for canonical OU gradients.
 *
 * <p>The likelihood delegate owns the canonical integrator lifecycle; this
 * adapter owns the gradient-specific parameter validation and dispatch.</p>
 */
public final class CanonicalOUGradientAdapter {

    private final CanonicalOUIntegrator integrator;
    private final int traitDimension;

    public CanonicalOUGradientAdapter(final CanonicalOUIntegrator integrator,
                                      final int traitDimension) {
        this.integrator = integrator;
        this.traitDimension = traitDimension;
    }

    public double[] getSelectionGradient(final Parameter requestedParameter,
                                         final AbstractBlockDiagonalTwoByTwoMatrixParameter nativeBlockParameter) {
        if (nativeBlockParameter == null) {
            throw new UnsupportedOperationException(
                    "Canonical OU XML wiring currently supports only native orthogonal-block selection gradients.");
        }
        if (requestedParameter.getDimension() != nativeBlockParameter.getParameter().getDimension()) {
            throw new IllegalArgumentException(
                    "Requested parameter dimension does not match native orthogonal-block parameter dimension.");
        }

        final double[] gradient = new double[requestedParameter.getDimension()];
        integrator.computeSelectionGradient(gradient);
        return gradient;
    }

    public double[] getMeanGradient(final Parameter requestedParameter,
                                    final boolean includeStationaryMean,
                                    final boolean includeRootMean) {
        if (!includeStationaryMean && !includeRootMean) {
            throw new IllegalArgumentException("Canonical mean gradient requested no active target.");
        }
        if (requestedParameter.getDimension() != traitDimension) {
            throw new IllegalArgumentException(
                    "Requested mean parameter dimension does not match canonical OU trait dimension.");
        }

        final double[] gradient = new double[requestedParameter.getDimension()];
        final double[] scratchGradient = new double[requestedParameter.getDimension()];

        if (includeStationaryMean) {
            integrator.computeStationaryMeanGradient(scratchGradient);
            addInPlace(gradient, scratchGradient);
            java.util.Arrays.fill(scratchGradient, 0.0);
        }
        if (includeRootMean) {
            integrator.computeGradientRootMean(scratchGradient);
            addInPlace(gradient, scratchGradient);
        }
        return gradient;
    }

    public double[] getPrecisionSourceGradient(final MatrixParameterInterface requestedParameter) {
        if (requestedParameter.getRowDimension() != traitDimension
                || requestedParameter.getColumnDimension() != traitDimension) {
            throw new IllegalArgumentException(
                    "Requested precision/variance matrix dimension does not match canonical OU trait dimension.");
        }

        final double[] gradient = new double[traitDimension * traitDimension];
        integrator.computeDiffusionGradient(gradient);
        return gradient;
    }

    private static void addInPlace(final double[] accumulator, final double[] increment) {
        for (int i = 0; i < accumulator.length; i++) {
            accumulator[i] += increment[i];
        }
    }
}
