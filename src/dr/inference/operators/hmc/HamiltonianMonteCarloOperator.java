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
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.math.distributions.NormalDistribution;
import dr.util.Transform;

/**
 * @author Max Tolkoff
 * @author Marc A. Suchard
 */

public class HamiltonianMonteCarloOperator extends AbstractCoercableOperator {

    final GradientWrtParameterProvider gradientProvider;
    protected double stepSize;
    protected final int nSteps;
    final NormalDistribution drawDistribution;
    final LeapFrogEngine leapFrogEngine;

    public HamiltonianMonteCarloOperator(CoercionMode mode, double weight, GradientWrtParameterProvider gradientProvider,
                                         Parameter parameter, Transform transform,
                                         double stepSize, int nSteps, double drawVariance) {
        super(mode);
        setWeight(weight);
        setTargetAcceptanceProbability(0.8); // Stan default

        this.gradientProvider = gradientProvider;
        this.stepSize = stepSize;
        this.nSteps = nSteps;
        this.drawDistribution = new NormalDistribution(0, Math.sqrt(drawVariance));
        this.leapFrogEngine = (transform != null ?
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

    static double getScaledDotProduct(final double[] momentum,
                                      final double sigmaSquared) {
        double total = 0.0;
        for (double m : momentum) {
            total += m * m;
        }

        return total / (2 * sigmaSquared);
    }

    static double[] drawInitialMomentum(final NormalDistribution distribution, final int dim) {
        double[] momentum = new double[dim];
        for (int i = 0; i < dim; i++) {
            momentum[i] = (Double) distribution.nextRandom();
        }
        return momentum;
    }

    @Override
    public double doOperation() { return leapFrog(); }

    private long count = 0;

    private static final boolean DEBUG = false;

    private double leapFrog() {

        if (DEBUG) {
            if (count % 5 == 0) {
                System.err.println(stepSize);
            }
            ++count;
        }

        final int dim = gradientProvider.getDimension();

        final double sigmaSquared = drawDistribution.getSD() * drawDistribution.getSD();

        final double[] momentum = drawInitialMomentum(drawDistribution, dim);
        final double[] position = leapFrogEngine.getInitialPosition();

        final double prop = getScaledDotProduct(momentum, sigmaSquared) +
                leapFrogEngine.getParameterLogJacobian();

        leapFrogEngine.updateMomentum(position, momentum,
                gradientProvider.getGradientLogDensity(), stepSize / 2);


        if (DEBUG) {
            System.err.println("nSteps = " + nSteps);
        }

        for (int i = 0; i < nSteps; i++) { // Leap-frog

            leapFrogEngine.updatePosition(position, momentum, stepSize, sigmaSquared);

            if (i < (nSteps - 1)) {
                leapFrogEngine.updateMomentum(position, momentum,
                        gradientProvider.getGradientLogDensity(), stepSize);
            }
        }

        leapFrogEngine.updateMomentum(position, momentum,
                gradientProvider.getGradientLogDensity(), stepSize / 2);

        final double res = getScaledDotProduct(momentum, sigmaSquared) +
                leapFrogEngine.getParameterLogJacobian();

        return prop - res; //hasting ratio
    }

    @Override
    public double getCoercableParameter() {
        return Math.log(stepSize);
    }

    @Override
    public void setCoercableParameter(double value) {
        if (DEBUG) {
            System.err.println("Setting coercable paramter: " + getCoercableParameter() + " -> " + value);
        }
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

        void setParameter(double[] position);

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
                }

                setParameter(position);
            }


            public void setParameter(double[] position) {

                final int dim = position.length;
                for (int j = 0; j < dim; ++j) {
                    parameter.setParameterValueQuietly(j, position[j]);
                }
                parameter.fireParameterChangedEvent();  // Does not seem to work with MaskedParameter
//                parameter.setParameterValueNotifyChangedAll(0, position[0]);
            }
        }

        class WithTransform extends Default {

            final private Transform transform;
            double[] unTransformedPosition;

            private WithTransform(Parameter parameter, Transform transform) {
                super(parameter);
                this.transform = transform;
            }

            @Override
            public double getParameterLogJacobian() {
                return transform.getLogJacobian(unTransformedPosition,0, unTransformedPosition.length);
            }

            @Override
            public double[] getInitialPosition() {
                unTransformedPosition = super.getInitialPosition();
                return transform.transform(unTransformedPosition, 0, unTransformedPosition.length);
            }

            @Override
            public void updateMomentum(double[] position, double[] momentum, double[] gradient,
                                       double functionalStepSize) {

                gradient = transform.updateGradientLogDensity(gradient, unTransformedPosition,
                        0, unTransformedPosition.length);

                super.updateMomentum(position, momentum, gradient, functionalStepSize);
            }

            @Override
            public void setParameter(double[] position) {
                unTransformedPosition = transform.inverse(position, 0, position.length);
                super.setParameter(unTransformedPosition);
            }
        }
    }
}
