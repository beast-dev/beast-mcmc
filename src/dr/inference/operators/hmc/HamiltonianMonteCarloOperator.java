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

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.GeneralOperator;
import dr.math.distributions.NormalDistribution;
import dr.util.Transform;

/**
 * @author Max Tolkoff
 * @author Marc A. Suchard
 */

public class HamiltonianMonteCarloOperator extends AbstractCoercableOperator {

    private final GradientWrtParameterProvider gradientProvider;
//    private final Parameter parameter;
//    private final Transform transform;
    private double stepSize;
    private final int nSteps;
    private final NormalDistribution drawDistribution;

    private final LeapFrogEngine leafFropEngine;

    public HamiltonianMonteCarloOperator(CoercionMode mode, double weight, GradientWrtParameterProvider gradientProvider,
                                         Parameter parameter, Transform transform,
                                         double stepSize, int nSteps, double drawVariance) {
        super(mode);
        setWeight(weight);
        setTargetAcceptanceProbability(0.8); // Stan default

        this.gradientProvider = gradientProvider;
//        this.parameter = parameter;
//        this.transform = transform;
        
        this.stepSize = stepSize;
        this.nSteps = nSteps;
        this.drawDistribution = new NormalDistribution(0, Math.sqrt(drawVariance));

        this.leafFropEngine = (transform != null ?
                new LeapFrogEngine.WithTransform(parameter, transform) :
                new LeapFrogEngine.Default(parameter));
    }

    @Override
    public String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public String getOperatorName() {
        return "Vanilla HMC operator";
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

    @Override
    public double doOperation() {
        return leafFrog();
    }

    protected double leafFrog() {

        final int dim = gradientProvider.getDimension();
        final double sigmaSquared = drawDistribution.getSD() * drawDistribution.getSD();

        double[] momentum = drawInitialMomentum(drawDistribution, dim);
        double[] position = leafFropEngine.getInitialPosition();

        final double prop = getScaledDotProduct(momentum, sigmaSquared) -
                leafFropEngine.getParameterLogJacobian();

        leafFropEngine.updateMomentum(position, momentum,
                gradientProvider.getGradientLogDensity(), stepSize / 2);

        for (int i = 0; i < nSteps; i++) { // Leap-frog

            leafFropEngine.updatePosition(position, momentum, stepSize, sigmaSquared);

            if (i < (nSteps - 1)) {
                leafFropEngine.updateMomentum(position, momentum,
                        gradientProvider.getGradientLogDensity(), stepSize);
            }
        } // end of loop over steps

        leafFropEngine.updateMomentum(position, momentum,
                gradientProvider.getGradientLogDensity(), stepSize / 2);

        final double res = getScaledDotProduct(momentum, sigmaSquared) -
                leafFropEngine.getParameterLogJacobian();

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

    interface LeapFrogEngine {

        double[] getInitialPosition();

        double getParameterLogJacobian();

        void updateMomentum(final double[] position,
                            final double[] momentum,
                            final double[] gradient,
                            final double functionalStepSize);

        void updatePosition(final double[] position,
                            final double[] momentum,
                            final double functionalStepSize,
                            final double sigmaSquared);

        class Default implements LeapFrogEngine {

            final protected Parameter parameter;

            protected Default(Parameter parameter) {
                this.parameter = parameter;
            }

            @Override
            public double[] getInitialPosition() {
                return parameter.getParameterValues();
            }

            @Override
            public double getParameterLogJacobian() {
                return 0;
            }

            @Override
            public void updateMomentum(double[] position, double[] momentum, double[] gradient,
                                       double functionalStepSize) {

                final int dim = momentum.length;
                for (int i = 0; i < dim; ++i) {
                    momentum[i] += functionalStepSize  * gradient[i];
                }
            }

            @Override
            public void updatePosition(double[] position, double[] momentum,
                                       double functionalStepSize, double sigmaSquared) {

                final int dim = momentum.length;
                for (int j = 0; j < dim; j++) {
                    position[j] += functionalStepSize * momentum[j] / sigmaSquared;

//                    final double oldValue = position[j];
//                    final double newValue = oldValue +
//                            functionalStepSize * momentum[j] / sigmaSquared;
//
//                    position[j] = newValue;
                }

                setParameter(position); // Write back into BEAST model
            }

            protected void setParameter(double[] position) {

                final int dim = position.length;
                for (int j = 0; j < dim; ++j) {
                    parameter.setParameterValueQuietly(j, position[j]);
                }
                parameter.fireParameterChangedEvent();
            }
        }

        class WithTransform extends Default {

            final private Transform transform;

            protected WithTransform(Parameter parameter, Transform transform) {
                super(parameter);
                this.transform = transform;
            }

            @Override
            public double getParameterLogJacobian() {
                return Math.log(parameter.getParameterValue(0)); // TODO
            }

            @Override
            public double[] getInitialPosition() {
                final double[] constrainedPosition = super.getInitialPosition();
                return transform.transform(constrainedPosition, 0, constrainedPosition.length);
            }

            @Override
            public void updateMomentum(double[] position, double[] momentum, double[] gradient, double functionalStepSize) {

                final int dim = momentum.length;


                    double[] x = parameter.getParameterValues(); // TODO
//            double[] transformGradient = transform.gradient(parameter.getParameterValues(), 0, dim);

                    for (int i = 0; i < dim; ++i) {
                        gradient[i] = gradient[i] * x[i] + 1.0; // TODO Trying to get log-transform to work
                    }

                    super.updateMomentum(position, momentum, gradient, functionalStepSize);
            }

            @Override
            protected void setParameter(double[] position) {
                double[] newConstrainedValues = transform.inverse(position, 0, position.length);
                super.setParameter(newConstrainedValues);
            }
        }
    }
}
