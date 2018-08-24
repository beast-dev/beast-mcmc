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

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.Vector;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Max Tolkoff
 * @author Marc A. Suchard
 */
public class SumDerivative implements GradientWrtParameterProvider {

    private final int dimension;
    private final Likelihood likelihood;
    private final Parameter parameter;

    List<GradientWrtParameterProvider> derivativeList;

    public SumDerivative(List<GradientWrtParameterProvider> derivativeList){

        // TODO Check that parameters are the same

        this.derivativeList = derivativeList;

        GradientWrtParameterProvider first = derivativeList.get(1);
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

        if (DEBUG) {
            // start timer
        }

        final double[] derivative = derivativeList.get(0).getGradientLogDensity();

        if (DEBUG) {
            // stop timer

            String name = derivativeList.get(0).getLikelihood().getId();
            System.err.println(name);
            System.err.println(new Vector(derivative));
        }

        for (int i = 1; i < size; i++) {

            if (DEBUG) {
                // start timer
            }

            final double[] temp = derivativeList.get(i).getGradientLogDensity();

            if (DEBUG) {
                // stop timer
                
                String name = derivativeList.get(i).getLikelihood().getId();
                System.err.println(name);
                System.err.println(new Vector(temp));
            }

            for (int j = 0; j < temp.length; j++) {
                derivative[j] += temp[j];
            }
        }

        if (DEBUG) {
            // print times

            if (DEBUG_KILL) {
                System.exit(-1);
            }
        }

        return derivative;
    }

    private static final boolean DEBUG = false;
    private static final boolean DEBUG_KILL = false;
}
