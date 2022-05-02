/*
 * CompoundGradient.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.inference.hmc;

import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */

@Deprecated
public class CompoundDerivative extends CompoundGradient implements HessianWrtParameterProvider{
    public CompoundDerivative(List<GradientWrtParameterProvider> derivativeList) {
        super(derivativeList);
    }

    @Override
    public double[] getDiagonalHessianLogDensity() {

        double[] result = new double[dimension];

        int offset = 0;
        for (GradientWrtParameterProvider derivative : derivativeList) {

            assert(derivative instanceof HessianWrtParameterProvider);

            double[] tmp = ((HessianWrtParameterProvider) derivative).getDiagonalHessianLogDensity();
            System.arraycopy(tmp, 0, result, offset, derivative.getDimension());
            offset += derivative.getDimension();
        }

        return result;
    }

    @Override
    public double[][] getHessianLogDensity() {
        throw new RuntimeException("Not implemented");
    }
}
