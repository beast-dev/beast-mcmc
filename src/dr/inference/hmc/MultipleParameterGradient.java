/*
 * MultipleParameterGradient.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;

import java.util.List;

/**
 * @author Alexander Fisher
 */

public class MultipleParameterGradient implements GradientWrtParameterProvider {
    //todo: include parameter order check with gradient order
    private final int dimension;
    private final Likelihood likelihood;
    private final Parameter parameter;
    private final int derivativeListSize;
    private final List<GradientWrtParameterProvider> derivativeList;
    private final int smallDim;

    public MultipleParameterGradient(List<GradientWrtParameterProvider> derivativeList, Parameter parameter) {
        this.derivativeList = derivativeList;
        this.derivativeListSize = derivativeList.size();
        this.dimension = parameter.getDimension();
        this.smallDim = dimension / derivativeListSize;

        // todo: 'likelihood' = prior, need to combine priors properly
        this.likelihood = derivativeList.get(0).getLikelihood();
        this.parameter = parameter;
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public int getDimension() {
        return dimension;
    }

    @Override
    public double[] getGradientLogDensity() {
        double[] fullDerivative = new double[dimension];
        for (int i = 0; i < derivativeListSize; i++) {
            System.arraycopy(getParameterGradientLogDensity(i), 0, fullDerivative, i * smallDim, smallDim);
        }
        return fullDerivative;
    }

    public double[] getParameterGradientLogDensity(int i) {
        return derivativeList.get(i).getGradientLogDensity();
    }
}
