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

import java.util.Arrays;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.math.MathUtils;
import dr.math.MultivariateFunction;
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
    private final double randomStepCountFraction;
    final NormalDistribution drawDistribution;
    final LeapFrogEngine leapFrogEngine;
    final MassProvider massProvider;

    public HamiltonianMonteCarloOperator(CoercionMode mode, double weight, GradientWrtParameterProvider gradientProvider,
                                         Parameter parameter, Transform transform,
                                         double stepSize, int nSteps, double drawVariance,
                                         double randomStepCountFraction) {
        this(mode, weight, gradientProvider, parameter, transform,
                stepSize, nSteps, drawVariance, randomStepCountFraction, false);
    }

    public HamiltonianMonteCarloOperator(CoercionMode mode, double weight, GradientWrtParameterProvider gradientProvider,
                                         Parameter parameter, Transform transform,
                                         double stepSize, int nSteps, double drawVariance,
                                         double randomStepCountFraction,
                                         boolean preConditioning) {
        super(mode);
        setWeight(weight);
        setTargetAcceptanceProbability(0.8); // Stan default

        this.gradientProvider = gradientProvider;
        this.stepSize = stepSize;
        this.nSteps = nSteps;
        this.randomStepCountFraction = randomStepCountFraction;
        this.drawDistribution = new NormalDistribution(0, 1.0);

        this.leapFrogEngine = (transform != null ?
                new LeapFrogEngine.WithTransform(parameter, transform, getDefaultInstabilityHandler()) :
                new LeapFrogEngine.Default(parameter, getDefaultInstabilityHandler()));

        this.massProvider = (!preConditioning ?
                new MassProvider.Default(gradientProvider.getDimension(), drawVariance) :
                (transform != null ?
                        new MassProvider.PreConditioningWithTransform(drawVariance, (HessianWrtParameterProvider) gradientProvider, transform) :
                        new MassProvider.PreConditioning(drawVariance, (HessianWrtParameterProvider) gradientProvider))
                );
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
                                      final double[] mass) {
        double total = 0.0;
        assert(momentum.length == mass.length);
        for (int i = 0; i < momentum.length; i++) {
            total += momentum[i] * momentum[i] / (mass[i] * 2.0);
        }
        return total;
    }

    static double[] drawInitialMomentum(final NormalDistribution distribution, final int dim, double[] mass) {
        double[] momentum = new double[dim];
        for (int i = 0; i < dim; i++) {
            momentum[i] = (Double) distribution.nextRandom() * Math.sqrt(mass[i]);
        }
        return momentum;
    }

    @Override
    public double doOperation() {
        try {
            return leapFrog();
        } catch (NumericInstabilityException e) {
            return Double.NEGATIVE_INFINITY;
        }
    }

    private long count = 0;

    private static final boolean DEBUG = false;
    
    static class NumericInstabilityException extends Exception { }

    private int getNumberOfSteps() {
        int count = nSteps;
        if (randomStepCountFraction > 0.0) {
            double draw = count * (1.0 + randomStepCountFraction * (MathUtils.nextDouble() - 0.5));
            count = Math.max(1, (int) draw);
        }
        return count;
    }

    private double leapFrog() throws NumericInstabilityException {

        if (DEBUG) {
            if (count % 5 == 0) {
                System.err.println("HMC step size: " + stepSize);
            }
            ++count;
        }

        final int dim = gradientProvider.getDimension();

        final double[] mass = massProvider.getMass();

        final double[] position = leapFrogEngine.getInitialPosition();
        final double[] momentum = drawInitialMomentum(drawDistribution, dim, mass);

        final double prop = getScaledDotProduct(momentum, mass) +
                leapFrogEngine.getParameterLogJacobian();

        leapFrogEngine.updateMomentum(position, momentum,
                gradientProvider.getGradientLogDensity(), stepSize / 2);

        int nStepsThisLeap = getNumberOfSteps();

        for (int i = 0; i < nStepsThisLeap; i++) { // Leap-frog

            leapFrogEngine.updatePosition(position, momentum, mass, stepSize);

            if (i < (nStepsThisLeap - 1)) {
                leapFrogEngine.updateMomentum(position, momentum,
                        gradientProvider.getGradientLogDensity(), stepSize);
            }
        }

        leapFrogEngine.updateMomentum(position, momentum,
                gradientProvider.getGradientLogDensity(), stepSize / 2);

        final double res = getScaledDotProduct(momentum, mass) +
                leapFrogEngine.getParameterLogJacobian();

        return prop - res; //hasting ratio
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

    enum InstabilityHandler {

        REJECT {
            @Override
            void checkValue(double x) throws NumericInstabilityException {
                if (Double.isNaN(x)) throw new NumericInstabilityException();
            }
        },

        DEBUG {
            @Override
            void checkValue(double x) throws NumericInstabilityException {
                if (Double.isNaN(x)) {
                    System.err.println("Numerical instability in HMC momentum; throwing exception");
                    throw new NumericInstabilityException();
                }
            }
        },

        IGNORE {
            @Override
            void checkValue(double x) {
                // Do nothing
            }
        };

        abstract void checkValue(double x) throws NumericInstabilityException;
    }

    protected InstabilityHandler getDefaultInstabilityHandler() {
        if (DEBUG) {
            return InstabilityHandler.DEBUG;
        } else {
            return InstabilityHandler.REJECT;
        }
    }

    interface LeapFrogEngine {

        double[] getInitialPosition();

        double getParameterLogJacobian();

        void updateMomentum(final double[] position,
                            final double[] momentum,
                            final double[] gradient,
                            final double functionalStepSize) throws NumericInstabilityException;

        void updatePosition(final double[] position,
                            final double[] momentum,
                            final double[] mass,
                            final double functionalStepSize);

        void setParameter(double[] position);

        class Default implements LeapFrogEngine {

            final protected Parameter parameter;
            final private InstabilityHandler instabilityHandler;

            protected Default(Parameter parameter, InstabilityHandler instabilityHandler) {
                this.parameter = parameter;
                this.instabilityHandler = instabilityHandler;
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
                                       double functionalStepSize) throws NumericInstabilityException {

                final int dim = momentum.length;
                for (int i = 0; i < dim; ++i) {
                    momentum[i] += functionalStepSize  * gradient[i];
                    instabilityHandler.checkValue(momentum[i]);
                }
            }

            @Override
            public void updatePosition(double[] position, double[] momentum, double[] mass,
                                       double functionalStepSize) {

                final int dim = momentum.length;
                for (int j = 0; j < dim; j++) {
                    position[j] += functionalStepSize * momentum[j] / mass[j];
                }

                setParameter(position);
            }

            public void setParameter(double[] position) {

                final int dim = position.length;
                for (int j = 0; j < dim; ++j) {
                    parameter.setParameterValueQuietly(j, position[j]);
                }
                parameter.fireParameterChangedEvent();  // Does not seem to work with MaskedParameter
            }
        }

        class WithTransform extends Default {

            final private Transform transform;
            double[] unTransformedPosition;

            private WithTransform(Parameter parameter, Transform transform, InstabilityHandler instabilityHandler) {
                super(parameter, instabilityHandler);
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
                                       double functionalStepSize) throws NumericInstabilityException {

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

    protected interface MassProvider {

        double[] getMass();

        class Default implements MassProvider {

            final int dim;
            final double[] mass;

            Default(int dim, double drawVariance) {
                this.dim = dim;
                this.mass = new double[dim];
                Arrays.fill(mass, drawVariance);
            }

            public double[] getMass() {
                return mass;
            }
        }

        class PreConditioning extends Default {

            final HessianWrtParameterProvider hessianWrtParameterProvider;

            PreConditioning(double drawVariance, HessianWrtParameterProvider hessianWrtParameterProvider) {
                super(hessianWrtParameterProvider.getDimension(), drawVariance);
                if (!(hessianWrtParameterProvider instanceof HessianWrtParameterProvider)) {
                    throw new IllegalArgumentException("Must provide a HessianProvider for preConditioning.");
                }
                this.hessianWrtParameterProvider = hessianWrtParameterProvider;
            }

            @Override
            public double[] getMass() {
                double[] diagonalHessian = hessianWrtParameterProvider.getDiagonalHessianLogDensity();
                return boundMass(diagonalHessian);
            }

            private double[] boundMass(double[] diagonalHessian) {

                double sum = 0.0;
                final double lowerBound = 1E-2;
                final double upperBound = 1E2;
                double[] boundedMass = new double[dim];

                for (int i = 0; i < dim; i++) {
                    boundedMass[i] = -diagonalHessian[i];
                    if (boundedMass[i] < lowerBound) {
                        boundedMass[i] = lowerBound;
                    } else if (boundedMass[i] > upperBound) {
                        boundedMass[i] = upperBound;
                    } else {
                        boundedMass[i] = -diagonalHessian[i];
                    }
                    sum += 1.0 / boundedMass[i];
                }
                final double mean = sum / dim;
                for (int i = 0; i < dim; i++) {
                    boundedMass[i] = mean * boundedMass[i];
                }
                return boundedMass;
            }
        }

        class PreConditioningWithTransform extends PreConditioning {

            final Transform transform;

            PreConditioningWithTransform(double drawVariance, HessianWrtParameterProvider hessianWrtParameterProvider,
                                         Transform transform) {
                super(drawVariance, hessianWrtParameterProvider);
                this.transform = transform;
            }

            @Override
            public double[] getMass() {
                double[] gradient = hessianWrtParameterProvider.getGradientLogDensity();
                double[] diagonalHessian = hessianWrtParameterProvider.getDiagonalHessianLogDensity();
                double[] unTransformedPosition = hessianWrtParameterProvider.getParameter().getParameterValues();
                diagonalHessian = transform.updateDiagonalHessianLogDensity(diagonalHessian, gradient, unTransformedPosition,
                        0, diagonalHessian.length);
//                double[] testHessian = NumericalDerivative.diagonalHessian(numeric1, transform.transform(hessianWrtParameterProvider.getParameter().getParameterValues(), 0, dim));
                return super.boundMass(diagonalHessian);
//                return setArbitraryMatrix(dim);
            }
//            private double[] setArbitraryMatrix(int dim) {
//                NormalDistribution drawNormal = new NormalDistribution(0.0, 1.0);
//                double multiplier = 1.0;
//                double[] mass = new double[dim];
//                for (int i = 0; i < dim; i++) {
//                    mass[i] = Math.abs((Double) drawNormal.nextRandom()) * multiplier;
//                    multiplier *= 1.1;
//                }
//                return mass;
//            }

//            private SymmetricMatrix getNumericalHessian() {
//                double[][] hessian = new double[dim][dim];
//                double[] oldUntransformedPosition = hessianWrtParameterProvider.getParameter().getParameterValues();
//                double[] oldTransformedPosition = transform.transform(oldUntransformedPosition, 0, dim);
//                double[][] gradientPlus = new double[dim][dim];
//                double[][] gradientMinus = new double[dim][dim];
//                double[] h = new double[dim];
//                for (int i = 0; i < dim; i++) {
//                    h[i] = MachineAccuracy.SQRT_SQRT_EPSILON*(Math.abs(oldTransformedPosition[i]) + 1.0);
//                    hessianWrtParameterProvider.getParameter().setParameterValue(i, Math.exp(oldTransformedPosition[i] + h[i]));
//                    double[] tempGradient = hessianWrtParameterProvider.getGradientLogDensity();
//                    gradientPlus[i] = transform.updateGradientLogDensity(tempGradient, hessianWrtParameterProvider.getParameter().getParameterValues(),
//                            0, dim);
//
//                    hessianWrtParameterProvider.getParameter().setParameterValue(i, Math.exp(oldTransformedPosition[i] - h[i]));
//                    tempGradient = hessianWrtParameterProvider.getGradientLogDensity();
//                    gradientMinus[i] = transform.updateGradientLogDensity(tempGradient, hessianWrtParameterProvider.getParameter().getParameterValues(),
//                            0, dim);
//                }
//                for (int i = 0; i < dim; i++) {
//                    for (int j = i; j < dim; j++) {
//                        hessian[j][i] = hessian[i][j] = (gradientPlus[i][j] - gradientMinus[i][j]) / (4.0 * h[j]) + (gradientPlus[j][i] - gradientMinus[j][i]) / (4.0 * h[i]);
//                    }
//                }
//                return new SymmetricMatrix(hessian);
//            }
//
            private MultivariateFunction numeric1 = new MultivariateFunction() {
                @Override
                public double evaluate(double[] argument) {

                    for (int i = 0; i < argument.length; ++i) {
                        hessianWrtParameterProvider.getParameter().setParameterValue(i, Math.exp(argument[i]));
                    }

//            treeDataLikelihood.makeDirty();
                    return hessianWrtParameterProvider.getLikelihood().getLogLikelihood();
                }

                @Override
                public int getNumArguments() {
                    return hessianWrtParameterProvider.getDimension();
                }

                @Override
                public double getLowerBound(int n) {
                    return Double.NEGATIVE_INFINITY;
                }

                @Override
                public double getUpperBound(int n) {
                    return Double.POSITIVE_INFINITY;
                }
            };
        }

    }
}
