/*
 * CanonicalDiagonalPrecisionGradient.java
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

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.adapter.CanonicalOUGradientAdapter;
import dr.inference.model.MatrixParameterInterface;

public final class CanonicalDiagonalPrecisionGradient extends DiagonalPrecisionGradient {

    private final CanonicalOUGradientAdapter gradientAdapter;
    private final MatrixParameterInterface sourceMatrixParameter;

    public CanonicalDiagonalPrecisionGradient(final TreeDataLikelihood likelihood,
                                              final ContinuousDataLikelihoodDelegate continuousData,
                                              final MatrixParameterInterface parameter) {
        super(new CanonicalGradientWrtPrecisionProvider(parameter.getColumnDimension()), likelihood, parameter);
        this.gradientAdapter = continuousData.getCanonicalOUGradientAdapter();
        this.sourceMatrixParameter = parameter;
    }

    @Override
    public double[] getGradientLogDensity() {
        return super.getGradientLogDensity(
                gradientAdapter.getPrecisionSourceGradient(sourceMatrixParameter));
    }
}
