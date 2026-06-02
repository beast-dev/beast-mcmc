/*
 * PathGradient.java
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

package dr.inference.hmc;

import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.PathDependent;

import java.util.Arrays;

/**
 * @author Marc A. Suchard
 */
public class PathGradient implements HessianWrtParameterProvider, PathDependent {

    private final int dimension;
    private final Likelihood likelihood;
    private final Parameter parameter;

    private final GradientWrtParameterProvider source;
    private final GradientWrtParameterProvider destination;

    private double beta = 1.0;

    public PathGradient(final GradientWrtParameterProvider source,
                        final GradientWrtParameterProvider destination){

        this.source = source;
        this.destination = destination;

        this.dimension = source.getDimension();
        this.parameter = source.getParameter();

//        if (destination.getParameter() != parameter) {
//            throw new RuntimeException("Invalid construction");
//        }
        if (destination.getDimension() != dimension) {
            throw new RuntimeException("Unequal parameter dimensions");
        }
        if (!Arrays.equals(destination.getParameter().getParameterValues(), parameter.getParameterValues())){
            throw new RuntimeException("Unequal parameter values");
        }

        this.likelihood = new Likelihood.Abstract(source.getLikelihood().getModel()) {

            @Override
            protected double calculateLogLikelihood() {

                double likelihood = source.getLikelihood().getLogLikelihood();

                if (beta != 1.0) {
                    likelihood = blend(likelihood, destination.getLikelihood().getLogLikelihood(), beta);
                }

                return likelihood;
            }
        };
    }

    @Override
    public void setPathParameter(double beta) {
        this.beta = beta;
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

        final double[] likelihood = source.getGradientLogDensity();

        if (beta != 1.0) {
            
            final double[] destination = this.destination.getGradientLogDensity();

            for (int i = 0; i < likelihood.length; ++i) {
                likelihood[i] = blend(likelihood[i], destination[i], beta);
            }
        }

        return likelihood;
    }

    private static double blend(double source, double destination, double beta) {
        return beta * source + (1.0 - beta) * destination;
    }

    @Override
    public double[] getDiagonalHessianLogDensity() {

        if (!(source instanceof HessianWrtParameterProvider) ||
                !(destination instanceof HessianWrtParameterProvider)) {
            throw new RuntimeException("Must use Hessian providers");
        }

        final double[] likelihood = ((HessianWrtParameterProvider) source).getDiagonalHessianLogDensity();

        if (beta != 1.0) {
            final double[] second = ((HessianWrtParameterProvider) destination).getDiagonalHessianLogDensity();

            for (int i = 0; i < likelihood.length; ++i) {
                likelihood[i] = blend(likelihood[i], second[i], beta);
            }
        }

        return likelihood;
    }

    @Override
    public double[][] getHessianLogDensity() {
        throw new RuntimeException("Not yet implemented");
    }
}
