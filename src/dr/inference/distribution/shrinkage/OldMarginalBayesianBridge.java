/*
 * OldMarginalBayesianBridge.java
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

package dr.inference.distribution.shrinkage;

import dr.inference.model.Parameter;
import dr.math.distributions.MarginalizedAlphaStableDistribution;

public class OldMarginalBayesianBridge extends OldBayesianBridgeLikelihood {

    public OldMarginalBayesianBridge(Parameter coefficients,
                                     Parameter globalScale,
                                     Parameter exponent) {
        super(coefficients, globalScale, exponent);
    }

    double calculateLogLikelihood() {

        final double scale = globalScale.getParameterValue(0);
        final double alpha = exponent.getParameterValue(0);

        double sum = 0.0;
        for (int i = 0; i < dim; ++i) {
            sum += MarginalizedAlphaStableDistribution.logPdf(coefficients.getParameterValue(i), scale, alpha);
        }
        return sum;
    }

    @Override
    double[] calculateGradientLogDensity() {

        final double scale = globalScale.getParameterValue(0);
        final double alpha = exponent.getParameterValue(0);

        double[] gradient = new double[dim];
        for (int i = 0; i < dim; ++i) {
            gradient[i] = MarginalizedAlphaStableDistribution.gradLogPdf(coefficients.getParameterValue(i),
                    scale, alpha);
        }
        return gradient;
    }

    @Override
    public Parameter getLocalScale() { return null; }
}