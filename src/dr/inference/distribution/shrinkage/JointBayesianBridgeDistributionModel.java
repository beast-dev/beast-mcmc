/*
 * JointBayesianBridgeDistributionModel.java
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
import dr.inference.model.PriorPreconditioningProvider;
import dr.math.distributions.NormalDistribution;

/**
 * @author Marc A. Suchard
 * @author Akihiko Nishimura
 */

public class JointBayesianBridgeDistributionModel extends BayesianBridgeDistributionModel
        implements PriorPreconditioningProvider {

    public JointBayesianBridgeDistributionModel(Parameter globalScale,
                                                Parameter localScale,
                                                Parameter exponent,
                                                Parameter slabWidth,
                                                int dim,
                                                boolean includeNormalizingConstant) {
        super(globalScale, exponent, dim, includeNormalizingConstant);
        this.localScale = localScale;
        this.slabWidth = slabWidth;

        if (dim != localScale.getDimension()) {
            throw new IllegalArgumentException("Invalid dimensions");
        }

        addVariable(localScale);
    }

    @Override
    public double getCoefficient(int i) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public Parameter getLocalScale() { return localScale; }

    @Override
    public Parameter getSlabWidth() { return slabWidth; }

    @Override
    double[] gradientLogPdf(double[] x) {

        double[] gradient = new double[dim];
        for (int i = 0; i < dim; ++i) {
            gradient[i] = NormalDistribution.gradLogPdf(x[i], 0, getStandardDeviation(i));
        }
        return gradient;
    }

    @Override
    public double logPdf(double[] x) {

        double pdf = 0.0;
        for (int i = 0; i < dim; ++i) {
            pdf += NormalDistribution.logPdf(x[i], 0, getStandardDeviation(i));
        }

        if (includeNormalizingConstant) {
            // TODO Add density of localScale variables
            throw new RuntimeException("Not yet implemented");
        }

        return pdf;
    }

    @Override
    public double getStandardDeviation(int index) {
        double globalLocalProduct = globalScale.getParameterValue(0) * localScale.getParameterValue(index);
        if (slabWidth != null) {
            double ratio = globalLocalProduct / slabWidth.getParameterValue(0);
            globalLocalProduct /= Math.sqrt(1.0 + ratio * ratio);
        }
        return globalLocalProduct;
    }

    @Override
    public double[] hessianLogPdf(double[] x) {

        double[] hessian = new double[dim];
        for (int i = 0; i < dim; ++i) {
            hessian[i] = NormalDistribution.hessianLogPdf(x[i], 0, getStandardDeviation(i));
        }
        return hessian;
    }

    @Override
    public double[] nextRandom() {
        double[] draws;
        if ( slabWidth != null) {
            draws = BayesianBridgeRNG.nextRandom(globalScale.getParameterValue(0), exponent.getParameterValue(0), slabWidth.getParameterValue(0), dim);
        } else {
            draws = BayesianBridgeRNG.nextRandom(globalScale.getParameterValue(0), exponent.getParameterValue(0), dim);
        }
        return draws;
    }

    private final Parameter localScale;
    private final Parameter slabWidth;
}