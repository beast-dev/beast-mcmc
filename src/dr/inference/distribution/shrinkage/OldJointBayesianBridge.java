/*
 * OldJointBayesianBridge.java
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
import dr.math.distributions.NormalDistribution;

public class OldJointBayesianBridge extends OldBayesianBridgeLikelihood {

    public OldJointBayesianBridge(Parameter coefficients,
                                  Parameter globalScale,
                                  Parameter localScale,
                                  Parameter exponent) {
        super(coefficients, globalScale, exponent);
        this.localScale = localScale;

        addVariable(localScale);
    }

    double calculateLogLikelihood() {
        double sum = 0.0;
        for (int i = 0; i < dim; ++i) {
            sum += NormalDistribution.logPdf(coefficients.getParameterValue(i), 0, getStandardDeviation(i));
        }

        // TODO Add density of localScale variables

        return sum;
    }

    @Override
    double[] calculateGradientLogDensity() {

        double[] gradient = new double[dim];
        for (int i = 0; i < dim; ++i) {
            gradient[i] = NormalDistribution.gradLogPdf(coefficients.getParameterValue(i),
                    0, getStandardDeviation(i));
        }
        return gradient;
    }

    @Override
    public Parameter getLocalScale() { return localScale; }

    private double getStandardDeviation(int index) {
        return globalScale.getParameterValue(0) * localScale.getParameterValue(index);
    }

    private final Parameter localScale;
}