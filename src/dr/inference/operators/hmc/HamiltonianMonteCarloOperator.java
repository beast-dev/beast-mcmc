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
import dr.math.MachineAccuracy;
import dr.math.MathUtils;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.SymmetricMatrix;
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
    final LeapFrogEngine leapFrogEngine;
    final MomentumProvider momentumProvider;

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

        this.leapFrogEngine = (transform != null ?
                new LeapFrogEngine.WithTransform(parameter, transform, getDefaultInstabilityHandler()) :
                new LeapFrogEngine.Default(parameter, getDefaultInstabilityHandler()));

        this.momentumProvider = (!preConditioning ?
                new MomentumProvider.Default(gradientProvider.getDimension(), drawVariance) :
                (transform != null ?
                        new MomentumProvider.PreConditioningWithTransform(drawVariance, (HessianWrtParameterProvider) gradientProvider, transform) :
                        new MomentumProvider.PreConditioning(drawVariance, (HessianWrtParameterProvider) gradientProvider))
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

        final double[] position = leapFrogEngine.getInitialPosition();
        final double[] momentum = momentumProvider.drawInitialMomentum();

        final double prop = momentumProvider.getScaledDotProduct(momentum) +
                leapFrogEngine.getParameterLogJacobian();

        leapFrogEngine.updateMomentum(position, momentum,
                gradientProvider.getGradientLogDensity(), stepSize / 2);

        int nStepsThisLeap = getNumberOfSteps();

        for (int i = 0; i < nStepsThisLeap; i++) { // Leap-frog

            leapFrogEngine.updatePosition(position, momentumProvider.weightMomentum(momentum), stepSize);

            if (i < (nStepsThisLeap - 1)) {
                leapFrogEngine.updateMomentum(position, momentum,
                        gradientProvider.getGradientLogDensity(), stepSize);
            }
        }

        leapFrogEngine.updateMomentum(position, momentum,
                gradientProvider.getGradientLogDensity(), stepSize / 2);

        final double res = momentumProvider.getScaledDotProduct(momentum) +
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
            public void updatePosition(double[] position, double[] momentum,
                                       double functionalStepSize) {

                final int dim = momentum.length;
                for (int j = 0; j < dim; j++) {
                    position[j] += functionalStepSize * momentum[j];
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

    protected interface MomentumProvider {

        @Deprecated
        double[] drawInitialMomentum();

        @Deprecated
        double getScaledDotProduct(double[] momentum);

        @Deprecated
        double[] weightMomentum(double[] momentum);

        double[] getMass();

        class Default implements MomentumProvider {
            final double drawVariance;
            final int dim;
            final double[] mass;
            MultivariateNormalDistribution drawDistribution;

            Default(int dim, double drawVariance) {
                this.dim = dim;
                this.drawVariance = drawVariance;
                this.drawDistribution = setDrawDistribution(drawVariance);

                this.mass = new double[dim];
                Arrays.fill(mass, drawVariance);
            }

            private MultivariateNormalDistribution setDrawDistribution(double drawVariance) {
                double[] mean = new double[dim];
                Arrays.fill(mean, 0.0);
                return new MultivariateNormalDistribution(mean, 1.0/drawVariance);
            }

            public double[] getMass() {
                return mass;
            }

            @Override
            public double[] drawInitialMomentum() {
                return (double[]) drawDistribution.nextRandom();
            }

            @Override
            public double getScaledDotProduct(double[] momentum) {
                double total = 0.0;
                for (int i = 0; i < momentum.length; i++) {
                    total += momentum[i] * momentum[i] / (2.0 * drawVariance);
                }
                return total;
            }

            @Override
            public double[] weightMomentum(double[] momentum) {
                double[] weightedMomentum = new double[dim];
                for (int i = 0; i < dim; i++) {
                    weightedMomentum[i] = momentum[i] / drawVariance;
                }
                return weightedMomentum;
            }
        }

        class PreConditioning extends Default {

            final HessianWrtParameterProvider hessianWrtParameterProvider;
            double[][] massMatrixInverse;

            PreConditioning(double drawVariance, HessianWrtParameterProvider hessianWrtParameterProvider) {
                super(hessianWrtParameterProvider.getDimension(), drawVariance);
                if (!(hessianWrtParameterProvider instanceof HessianWrtParameterProvider)) {
                    throw new IllegalArgumentException("Must provide a HessianProvider for preConditioning.");
                }
                this.hessianWrtParameterProvider = hessianWrtParameterProvider;
                this.massMatrixInverse = new double[dim][dim];
//                setMassMatrixInverse(hessianWrtParameterProvider.getDiagonalHessianLogDensity());
            }

            public void setMassMatrixInverse(double[] diagonalHessian) {
                for (int i = 0; i < dim; i++) {
                    Arrays.fill(massMatrixInverse[i], 0.0);
                }
                boundSigmaSquaredInverse(diagonalHessian, drawVariance);
                for (int i = 0; i < dim; i++) {
                    massMatrixInverse[i][i] = diagonalHessian[i];
                }
            }

            private MultivariateNormalDistribution setDrawDistribution(double drawVariance) {
                double[] mean = new double[dim];
                Arrays.fill(mean, 0.0);
                return new MultivariateNormalDistribution(mean, massMatrixInverse);
            }

            @Override
            public double[] drawInitialMomentum() {
                updateMassMatrixInverse();
                return (double[]) drawDistribution.nextRandom();
            }

            private void updateMassMatrixInverse() {
                setMassMatrixInverse(hessianWrtParameterProvider.getDiagonalHessianLogDensity());
                this.drawDistribution = setDrawDistribution(drawVariance);
            }

            @Override
            public double getScaledDotProduct(double[] momentum) {
                double total = 0.0;
                for (int i = 0; i < momentum.length; i++) {
                    total += momentum[i] * momentum[i] * massMatrixInverse[i][i] / 2.0;
                }
                return total;
            }

            @Override
            public double[] weightMomentum(double[] momentum) {
                double[] weightedMomentum = new double[dim];
                for (int i = 0; i < dim; i++) {
                    weightedMomentum[i] = massMatrixInverse[i][i] * momentum[i];
                }
                return weightedMomentum;
            }

//            private void boundMass(double[] mass) {
//
//            }
//
//            public double[] getMass() {
//                return null;
//            }

            private void boundSigmaSquaredInverse(double[] sigmaSquaredInverse, double drawVariance) {

                double min, max, mean, sum;
                min = max =  Math.log(Math.abs(sigmaSquaredInverse[0]) / drawVariance);
                for (int i = 0; i < sigmaSquaredInverse.length; i++) {
                    sigmaSquaredInverse[i] = Math.log(Math.abs(sigmaSquaredInverse[i]) / drawVariance);
                    if (sigmaSquaredInverse[i] > max) max = sigmaSquaredInverse[i];
                    if (sigmaSquaredInverse[i] < min) min = sigmaSquaredInverse[i];
                }
                sum = 0.0;
                if (max - min > 4.0) {
                    for (int i = 0; i < sigmaSquaredInverse.length; i++) {
                        sigmaSquaredInverse[i] = Math.exp(-((sigmaSquaredInverse[i] - min) / (max - min) * 4.0 - 2.0));
                        sum += sigmaSquaredInverse[i];
                    }
                } else {
                    for (int i = 0; i < sigmaSquaredInverse.length; i++) {
                        sigmaSquaredInverse[i] = Math.exp(-sigmaSquaredInverse[i]); //Math.exp(-((sigmaSquaredInverse[i] - min) / (max - min) * 4.0 - 2.0));
                        sum += sigmaSquaredInverse[i];
                    }
                }

                mean = sum / sigmaSquaredInverse.length;
                for (int i = 0; i < sigmaSquaredInverse.length; i++) {
                    sigmaSquaredInverse[i] /= mean;
                }
            }
        }

        class PreConditioningWithTransform extends PreConditioning {

            final Transform transform;

            PreConditioningWithTransform(double drawVariance, HessianWrtParameterProvider hessianWrtParameterProvider,
                                         Transform transform) {
                super(drawVariance, hessianWrtParameterProvider);
                this.transform = transform;
            }

            private void setArbitraryMatrix(double[][] matrix) {
                double[] diagonal = drawDistribution.nextMultivariateNormal();
                double multiplier = 1.2;
                for (int i = 0; i < dim; i++) {
                    matrix[i][i] = Math.abs(diagonal[i]) * multiplier;
                }
            }

            public void setMassMatrixInverse(double[] diagonalHessian) {
//                double[] gradient = hessianWrtParameterProvider.getGradientLogDensity();
//                double[] unTransformedPosition = hessianWrtParameterProvider.getParameter().getParameterValues();
//                diagonalHessian = transform.updateDiagonalHessianLogDensity(diagonalHessian, gradient, unTransformedPosition,
//                        0, diagonalHessian.length);
//
//                double[] testHessian = NumericalDerivative.diagonalHessian(numeric1, transform.transform(hessianWrtParameterProvider.getParameter().getParameterValues(), 0, dim));
                SymmetricMatrix hessian = getNumericalHessian();
                double[][] hessianInverse = hessian.inverse().toComponents();


//                super.setMassMatrixInverse(diagonalHessian);
//                setArbitraryMatrix(massMatrixInverse);
                for (int i = 0; i < dim; i++) {
                    for (int j = 0; j < dim; j++) {
                        massMatrixInverse[i][j] = -hessianInverse[i][j];
                    }
                }
                double[] diagonalInverseHessian = new double[dim];
                for (int i = 0; i < dim; i++) {
                    diagonalInverseHessian[i] = hessianInverse[i][i];
                }
                super.boundSigmaSquaredInverse(diagonalInverseHessian, drawVariance);
                for (int i = 0; i < dim; i++) {
                    massMatrixInverse[i][i] = diagonalInverseHessian[i];
                }
            }

            @Override
            public double[] weightMomentum(double[] momentum) {
                double[] weightedMomentum = new double[dim];
                for (int i = 0; i < dim; i++) {
                    double sum = 0;
                    for (int j = 0; j < dim; j++) {
                        sum += massMatrixInverse[i][j] * momentum[j];
                    }
                    weightedMomentum[i] = sum;
                }
                return weightedMomentum;
            }

            private SymmetricMatrix getNumericalHessian() {
                double[][] hessian = new double[dim][dim];
                double[] oldUntransformedPosition = hessianWrtParameterProvider.getParameter().getParameterValues();
                double[] oldTransformedPosition = transform.transform(oldUntransformedPosition, 0, dim);
                double[][] gradientPlus = new double[dim][dim];
                double[][] gradientMinus = new double[dim][dim];
                double[] h = new double[dim];
                for (int i = 0; i < dim; i++) {
                    h[i] = MachineAccuracy.SQRT_SQRT_EPSILON*(Math.abs(oldTransformedPosition[i]) + 1.0);
                    hessianWrtParameterProvider.getParameter().setParameterValue(i, Math.exp(oldTransformedPosition[i] + h[i]));
                    double[] tempGradient = hessianWrtParameterProvider.getGradientLogDensity();
                    gradientPlus[i] = transform.updateGradientLogDensity(tempGradient, hessianWrtParameterProvider.getParameter().getParameterValues(),
                            0, dim);

                    hessianWrtParameterProvider.getParameter().setParameterValue(i, Math.exp(oldTransformedPosition[i] - h[i]));
                    tempGradient = hessianWrtParameterProvider.getGradientLogDensity();
                    gradientMinus[i] = transform.updateGradientLogDensity(tempGradient, hessianWrtParameterProvider.getParameter().getParameterValues(),
                            0, dim);
                }
                for (int i = 0; i < dim; i++) {
                    for (int j = i; j < dim; j++) {
                        hessian[j][i] = hessian[i][j] = (gradientPlus[i][j] - gradientMinus[i][j]) / (4.0 * h[j]) + (gradientPlus[j][i] - gradientMinus[j][i]) / (4.0 * h[i]);
                    }
                }
                return new SymmetricMatrix(hessian);
            }

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
