/*
 * NewHamiltonianMonteCarloOperator.java
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

package dr.inference.operators.hmc;

import dr.inference.model.GradientWrtParameterProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.math.distributions.NormalDistribution;

/**
 * @author Max Tolkoff
 * @author Marc A. Suchard
 */

public class NewHamiltonianMonteCarloOperator extends AbstractCoercableOperator {

    private final GradientWrtParameterProvider gradientProvider;
    private final Parameter parameter;
    private double stepSize;
    private final int nSteps;
    private final NormalDistribution drawDistribution;

    public NewHamiltonianMonteCarloOperator(CoercionMode mode, double weight, GradientWrtParameterProvider gradientProvider,
                                            Parameter parameter, double stepSize, int nSteps, double drawVariance) {
        super(mode);
        setWeight(weight);

        this.gradientProvider = gradientProvider;
        this.parameter = parameter;
        this.stepSize = stepSize;
        this.nSteps = nSteps;
        this.drawDistribution = new NormalDistribution(0, Math.sqrt(drawVariance));
    }

    @Override
    public String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public String getOperatorName() {
        return "Vanilla HMC operator";
    }

    private static void updateHalfMomentum(final double[] momentum,
                                           final double functionalStepSize,
                                           final double[] gradient) {
        final int dim = momentum.length;

        for (int i = 0; i < dim; ++i) {
            momentum[i] = momentum[i] + functionalStepSize / 2 * gradient[i]; /* Sign change */
        }
    }

    private static double getScaledDotProduct(final double[] momentum,
                                              final double sigmaSquared) {
        final int dim = momentum.length;

        double total = 0.0;
        for (int i = 0; i < dim; i++) {
            total += momentum[i] * momentum[i];
        }

        return total / (2 * sigmaSquared);
    }

    private static double[] drawInitialMomentum(final NormalDistribution distribution, final int dim) {
        double[] momentum = new double[dim];
        for (int i = 0; i < dim; i++) {
            momentum[i] = (Double) distribution.nextRandom();
        }
        return momentum;
    }

    private static void updateParameter(final Parameter parameter,
                                        final double[] momentum,
                                        final double functionalStepSize,
                                        final double sigmaSquared) {
        final int dim = momentum.length;

        for (int j = 0; j < dim; j++) {
            final double oldValue = parameter.getParameterValue(j);
            final double newValue =  oldValue +
                    functionalStepSize * momentum[j] / sigmaSquared;

            parameter.setParameterValueQuietly(j, newValue);
        }
        parameter.fireParameterChangedEvent();
    }

    @Override
    public double doOperation() {

        final double functionalStepSize = stepSize;
        final int dim = gradientProvider.getDimension();
        final double sigmaSquared = drawDistribution.getSD() * drawDistribution.getSD();

        double[] momentum = drawInitialMomentum(drawDistribution, dim);
        final double prop = getScaledDotProduct(momentum, sigmaSquared);

        for (int i = 0; i < nSteps; i++) { // Leap-frog

            updateHalfMomentum(momentum, functionalStepSize,
                    gradientProvider.getGradientLogDensity());

            updateParameter(parameter, momentum, functionalStepSize, sigmaSquared);

            updateHalfMomentum(momentum, functionalStepSize,
                    gradientProvider.getGradientLogDensity());

        } // end of loop over steps

        final double res = getScaledDotProduct(momentum, sigmaSquared);

        return prop - res;
    }

    @Override
    public double getCoercableParameter() {
        return Math.log(stepSize);
    }

    @Override
    public void setCoercableParameter(double value) {
        stepSize = Math.exp(value);
    }

    @Override
    public double getRawParameter() {
        return stepSize;
    }
}
