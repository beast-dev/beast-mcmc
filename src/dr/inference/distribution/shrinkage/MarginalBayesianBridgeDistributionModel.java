/*
 * MarginalBayesianBridgeDistributionModel.java
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
import dr.math.distributions.LaplaceDistribution;
import dr.math.distributions.MarginalizedAlphaStableDistribution;

/**
 * @author Marc A. Suchard
 * @author Akihiko Nishimura
 */

public class MarginalBayesianBridgeDistributionModel extends BayesianBridgeDistributionModel {

    public MarginalBayesianBridgeDistributionModel(Parameter globalScale,
                                                   Parameter exponent,
                                                   int dim,
                                                   boolean includeNormalizingConstant) {
        super(globalScale, exponent, dim, includeNormalizingConstant);
    }

    @Override
    public double getCoefficient(int i) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public Parameter getLocalScale() {
        return null;
    }

    @Override
    public Parameter getSlabWidth() {
        return null;
    }

    @Override
    double[] gradientLogPdf(double[] x) {
        final int dim = x.length;
        final double scale = globalScale.getParameterValue(0);
        final double alpha = exponent.getParameterValue(0);

        double[] gradient = new double[dim];
        if (alpha != 1.0) {
            for (int i = 0; i < dim; ++i) {
                gradient[i] = MarginalizedAlphaStableDistribution.gradLogPdf(x[i], scale, alpha);
            }
        } else {
            for (int i = 0; i < dim; ++i) {
                gradient[i] = LaplaceDistribution.gradLogPdf(x[i], 0, scale);
            }
        }

        return gradient;
    }

    @Override
    public double logPdf(double[] v) {
        final double scale = globalScale.getParameterValue(0);
        final double alpha = exponent.getParameterValue(0);

        double sum = 0.0;
        if (alpha != 1.0) {
            for (double x : v) {
                sum += MarginalizedAlphaStableDistribution.logPdf(x, scale, alpha);
            }
        } else {
            for (int i = 0; i < dim; ++i) {
                sum += LaplaceDistribution.logPdf(v[i], 0, scale);
            }
        }

        if (includeNormalizingConstant) {
            // TODO Add
            throw new RuntimeException("Not yet implemented");
        }

        return sum;
    }

    @Override
    public double[] hessianLogPdf(double[] x) {
        throw new RuntimeException("Not yet implemented");
    }
}