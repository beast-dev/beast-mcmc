/*
 * SumDerivative.java
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

import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Max Tolkoff
 * @author Marc A. Suchard
 */
public class SumDerivative implements GradientWrtParameterProvider {

    private final int dimension;
    private final Likelihood likelihood;
    private final Parameter parameter;

    private final List<GradientWrtParameterProvider> derivativeList;

    public SumDerivative(List<GradientWrtParameterProvider> derivativeList){

        this.derivativeList = derivativeList;

        GradientWrtParameterProvider first = derivativeList.get(0);
        dimension = first.getDimension();
        parameter = first.getParameter();

        if (derivativeList.size() == 1) {
            likelihood = first.getLikelihood();
        } else {
            List<Likelihood> likelihoodList = new ArrayList<Likelihood>();

            for (GradientWrtParameterProvider grad : derivativeList) {
                if (grad.getDimension() != dimension) {
                    throw new RuntimeException("Unequal parameter dimensions");
                }
                if (!Arrays.equals(grad.getParameter().getParameterValues(), parameter.getParameterValues())){
                    throw new RuntimeException("Unequal parameter values");
                }
                likelihoodList.add(grad.getLikelihood());
            }
            likelihood = new CompoundLikelihood(likelihoodList);
        }
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
        int size = derivativeList.size();

        final double[] derivative = derivativeList.get(0).getGradientLogDensity();

        for (int i = 1; i < size; i++) {
            
            final double[] temp = derivativeList.get(i).getGradientLogDensity();

            for (int j = 0; j < temp.length; j++) {
                derivative[j] += temp[j];
            }
        }

        return derivative;
    }
}
