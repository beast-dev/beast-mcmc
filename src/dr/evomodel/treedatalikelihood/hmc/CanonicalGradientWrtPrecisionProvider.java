/*
 * CanonicalGradientWrtPrecisionProvider.java
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

package dr.evomodel.treedatalikelihood.hmc;

/**
 * Canonical OU precision-gradient provider backed by an externally supplied
 * covariance-space source gradient.
 *
 * <p>The canonical tree passer computes a flattened gradient with respect to the
 * diffusion covariance matrix {@code Q}. This provider reuses the existing
 * precision/variance chain rule from the HMC layer without going through the
 * legacy branch-specific gradient path.
 */
public final class CanonicalGradientWrtPrecisionProvider
        extends GradientWrtPrecisionProvider.AbstractGradientWrtPrecisionProvider {

    private final double[] symmetricGradientScratch;

    public CanonicalGradientWrtPrecisionProvider(final int dim) {
        this.dim = dim;
        this.symmetricGradientScratch = new double[dim * dim];
    }

    @Override
    public double[] getGradientWrtPrecision(final double[] vecV, final double[] gradient) {
        symmetrizeInto(gradient, symmetricGradientScratch);
        final MultivariateChainRule rule = new MultivariateChainRule.InverseGeneral(vecV);
        return rule.chainGradient(symmetricGradientScratch);
    }

    @Override
    public double[] getGradientWrtVariance(final double[] vecP, final double[] vecV, final double[] gradient) {
        symmetrizeInto(gradient, symmetricGradientScratch);
        return symmetricGradientScratch.clone();
    }

    private void symmetrizeInto(final double[] source, final double[] out) {
        for (int row = 0; row < dim; row++) {
            final int rowOffset = row * dim;
            for (int col = 0; col < dim; col++) {
                out[rowOffset + col] = 0.5 * (source[rowOffset + col] + source[col * dim + row]);
            }
        }
    }
}
