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
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.math.distributions.NormalDistribution;
import dr.util.Transform;

/**
 * @author Max Tolkoff
 * @author Marc A. Suchard
 */

public class HamiltonianMonteCarloOperator extends AbstractCoercableOperator {

    private final GradientWrtParameterProvider gradientProvider;
    private double stepSize;
    private final int nSteps;
    private final NormalDistribution drawDistribution;
    private final int[] mask;

    private final LeapFrogEngine leafFropEngine;

    public HamiltonianMonteCarloOperator(CoercionMode mode, double weight, GradientWrtParameterProvider gradientProvider,
                                         Parameter parameter, Transform transform, Parameter maskParameter,
                                         double stepSize, int nSteps, double drawVariance) {
        super(mode);
        setWeight(weight);
        setTargetAcceptanceProbability(0.8); // Stan default

        this.gradientProvider = gradientProvider;
        this.stepSize = stepSize;
        this.nSteps = nSteps;
        this.drawDistribution = new NormalDistribution(0, Math.sqrt(drawVariance));

        if (maskParameter != null) {
            mask = new int[maskParameter.getDimension()];
            for (int i = 0; i < maskParameter.getDimension(); ++i) {
                if (maskParameter.getParameterValue(i) == 1.0) {
                    mask[i] = 1;
                }
            }
        } else {
            mask = null;
        }

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

    private static double[] drawInitialMomentum(final NormalDistribution distribution, final int dim,
                                                final int[] mask) {
        double[] momentum = new double[dim];
        for (int i = 0; i < dim; i++) {
            if (mask == null || mask[i] == 1) {
                momentum[i] = (Double) distribution.nextRandom();
            }
        }
        return momentum;
    }

    @Override
    public double doOperation() {
        return leafFrog();
    }

    private long count = 0;

    protected double leafFrog() {

        if (count % 10000 == 0) {
            System.err.println(stepSize);
        }
        ++count;

        final int dim = gradientProvider.getDimension();
        final double sigmaSquared = drawDistribution.getSD() * drawDistribution.getSD();

        final double[] momentum = drawInitialMomentum(drawDistribution, dim, mask);
        final double[] position = leafFropEngine.getInitialPosition();

        final double prop = getScaledDotProduct(momentum, sigmaSquared) +
                leafFropEngine.getParameterLogJacobian();

        leafFropEngine.updateMomentum(position, momentum,
                gradientProvider.getGradientLogDensity(), stepSize / 2, mask);

        for (int i = 0; i < nSteps; i++) { // Leap-frog

            leafFropEngine.updatePosition(position, momentum, stepSize, sigmaSquared);

            if (i < (nSteps - 1)) {
                leafFropEngine.updateMomentum(position, momentum,
                        gradientProvider.getGradientLogDensity(), stepSize, mask);
            }
        } // end of loop over steps

        leafFropEngine.updateMomentum(position, momentum,
                gradientProvider.getGradientLogDensity(), stepSize / 2, mask);

        final double res = getScaledDotProduct(momentum, sigmaSquared) +
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
                            final double functionalStepSize,
                            final int[] mask);

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
                                       double functionalStepSize, int[] mask) {
                
                if (mask != null) {
                    for (int i = 0; i < mask.length; ++i) {
                        if (mask[i] == 0) {
                            gradient[i] = 0.0;
                        }
                    }
                }

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
            double[] unTransformedPosition;

            protected WithTransform(Parameter parameter, Transform transform) {
                super(parameter);
                this.transform = transform;
            }

            @Override
            public double getParameterLogJacobian() {
                return transform.getLogJacobian(unTransformedPosition, 0, unTransformedPosition.length);
            }

            @Override
            public double[] getInitialPosition() {
                unTransformedPosition = super.getInitialPosition();
                return transform.transform(unTransformedPosition, 0, unTransformedPosition.length);
            }

            @Override
            public void updateMomentum(double[] position, double[] momentum, double[] gradient,
                                       double functionalStepSize, int[] mask) {

                gradient = transform.updateGradientLogDensity(gradient, unTransformedPosition,
                        0, unTransformedPosition.length);

                super.updateMomentum(position, momentum, gradient, functionalStepSize, mask);
            }

            @Override
            protected void setParameter(double[] position) {
                unTransformedPosition = transform.inverse(position, 0, position.length);
                super.setParameter(unTransformedPosition);
            }
        }
    }
}
