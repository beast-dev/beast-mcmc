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

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Alexander Fisher
 */

public class MultipleParameterGradient implements GradientWrtParameterProvider {
    //todo: include parameter order check with gradient order
    private final int dimension;
    private final Likelihood likelihood;
    private final Parameter parameter;
    private final List<GradientWrtParameterProvider> derivativeList;

    public MultipleParameterGradient(List<GradientWrtParameterProvider> derivativeList, Parameter parameter) {
        this.derivativeList = derivativeList;
        int listSize = derivativeList.size();
        int totalDim = 0;
        // todo: remove since it's redundant
        for (int i = 0; i < listSize; i++) {
            totalDim = totalDim + derivativeList.get(i).getDimension();
        }
        if (totalDim != parameter.getDimension()) {
            throw new RuntimeException("Parameter dimension mismatch");
        }
        this.dimension = totalDim;
        // todo: check for same likelihood across derivativeList in parser
        this.likelihood = null;
        this.parameter = parameter;

//        this.derivativeList = derivativeList;
//
//        GradientWrtParameterProvider first = derivativeList.get(0);
//        dimension = first.getDimension();
//        parameter = first.getParameter();
//
//        if (derivativeList.size() == 1) {
//            likelihood = first.getLikelihood();
//        } else {
//            List<Likelihood> likelihoodList = new ArrayList<>();
//
//            for (GradientWrtParameterProvider grad : derivativeList) {
//                if (grad.getDimension() != dimension) {
//                    throw new RuntimeException("Unequal parameter dimensions");
//                }
//                if (!Arrays.equals(grad.getParameter().getParameterValues(), parameter.getParameterValues())){
//                    throw new RuntimeException("Unequal parameter values");
//                }
//                for (Likelihood likelihood : grad.getLikelihood().getLikelihoodSet()) {
//                    if (!(likelihoodList.contains(likelihood))) {
//                        likelihoodList.add(likelihood);
//                    }
//                }
//            }
//            likelihood = new CompoundLikelihood(likelihoodList);
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
        return new double[0];
    }

    public double getParameterGradientLogDensity(int i) {
        return 2.0;
    }
}
