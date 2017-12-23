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

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Max Tolkoff
 * @author Marc A. Suchard
 */

public class CompoundGradient implements GradientWrtParameterProvider {

    private final int dimension;
    private final List<GradientWrtParameterProvider> derivativeList;
    private final Likelihood likelihood;
    private final Parameter parameter;

    public CompoundGradient(List<GradientWrtParameterProvider> derivativeList) {

        this.derivativeList = derivativeList;

        if (derivativeList.size() == 1) {
            likelihood = derivativeList.get(0).getLikelihood();
            parameter = derivativeList.get(0).getParameter();
            dimension = parameter.getDimension();
        } else {
            List<Likelihood> likelihoodList = new ArrayList<Likelihood>();
            CompoundParameter compoundParameter = new CompoundParameter("hmc");

            int dim = 0;
            for (GradientWrtParameterProvider grad : derivativeList) {
                likelihoodList.add(grad.getLikelihood());
                Parameter p = grad.getParameter();
                compoundParameter.addParameter(p);

                dim += p.getDimension();
            }

            likelihood = new CompoundLikelihood(likelihoodList);
            parameter = compoundParameter;
            dimension = dim;
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

//    @Override
//    public void getGradientLogDensity(final double[] destination, final int offset) {
//        double[] grad = getGradientLogDensity();
//        System.arraycopy(grad, 0, destination, offset, grad.length);
//    }

    @Override
    public double[] getGradientLogDensity() {

        double[] result = new double[dimension];

        int offset = 0;
        for (GradientWrtParameterProvider grad : derivativeList) {

            System.err.println(grad.getLikelihood().getId() + " " + grad.getParameter().getId());

            double[] tmp = grad.getGradientLogDensity();
            System.arraycopy(tmp, 0, result, offset, grad.getDimension());
            offset += grad.getDimension();
        }

        return result;
    }
}
