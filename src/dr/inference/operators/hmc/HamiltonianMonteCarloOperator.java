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
    final boolean preConditioning;
    final double drawVariance;
    double[] sigmaSquaredInverse;
    double[] sigmaList;

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
        this.preConditioning = preConditioning;
        if (preConditioning && !(gradientProvider instanceof HessianWrtParameterProvider)) {
            throw new IllegalArgumentException("Must provide a HessianProvider for preConditioning.");
        }
        sigmaSquaredInverse = new double[gradientProvider.getDimension()];
        sigmaList = new double[gradientProvider.getDimension()];
        drawDistribution = new NormalDistribution(0, 1.0);

        this.drawVariance = drawVariance;
        setSigmas();
    }

    private void setSigmas() {
        if (preConditioning) {
            sigmaSquaredInverse = ((HessianWrtParameterProvider) gradientProvider).getDiagonalHessianLogDensity();
            boundSigmaInverse(sigmaSquaredInverse);
            for (int i = 0; i < sigmaSquaredInverse.length; i++) {
                sigmaList[i] = Math.sqrt(1.0 / sigmaSquaredInverse[i]);
            }
        } else {
            Arrays.fill(sigmaSquaredInverse, 1.0 / drawVariance);
            Arrays.fill(sigmaList, Math.sqrt(drawVariance));
        }
    }

    private void boundSigmaInverse(double[] sigmaSquaredInverse) {

        double min, max, mean, sum;
        min = max =  Math.log(Math.abs(sigmaSquaredInverse[0]) / drawVariance);
        for (int i = 0; i < sigmaSquaredInverse.length; i++) {
            sigmaSquaredInverse[i] = Math.log(Math.abs(sigmaSquaredInverse[i]) / drawVariance);
            if (sigmaSquaredInverse[i] > max) max = sigmaSquaredInverse[i];
            if (sigmaSquaredInverse[i] < min) min = sigmaSquaredInverse[i];
        }
        sum = 0.0;
        for (int i = 0; i < sigmaSquaredInverse.length; i++) {
            sigmaSquaredInverse[i] = Math.exp((sigmaSquaredInverse[i] - min) / (max - min) * 4.0 - 2.0);
            sum += sigmaSquaredInverse[i];
        }
        mean = sum / sigmaSquaredInverse.length;
        for (int i = 0; i < sigmaSquaredInverse.length; i++) {
            sigmaSquaredInverse[i] /= mean;
        }
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
                                      final double[] sigmaSquaredInverse) {
        double total = 0.0;
        assert(momentum.length == sigmaSquaredInverse.length);
        for (int i = 0; i < momentum.length; i++) {
            total += momentum[i] * momentum[i] * sigmaSquaredInverse[i] / 2.0;
        }

        return total;
    }

    static double[] drawInitialMomentum(final NormalDistribution distribution, final int dim, final double[] sigmaList) {
        double[] momentum = new double[dim];
        for (int i = 0; i < dim; i++) {
            momentum[i] = (Double) distribution.nextRandom() * sigmaList[i];
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

//        final double sigmaSquared = drawDistribution.getSD() * drawDistribution.getSD();
        if (preConditioning) {
            setSigmas();
        }

        final double[] momentum = drawInitialMomentum(drawDistribution, dim, sigmaList);
        final double[] position = leapFrogEngine.getInitialPosition();

        final double prop = getScaledDotProduct(momentum, sigmaSquaredInverse) +
                leapFrogEngine.getParameterLogJacobian();

        leapFrogEngine.updateMomentum(position, momentum,
                gradientProvider.getGradientLogDensity(), stepSize / 2);

        int nStepsThisLeap = getNumberOfSteps();

        for (int i = 0; i < nStepsThisLeap; i++) { // Leap-frog

            leapFrogEngine.updatePosition(position, momentum, stepSize, sigmaSquaredInverse);

            if (i < (nStepsThisLeap - 1)) {
                leapFrogEngine.updateMomentum(position, momentum,
                        gradientProvider.getGradientLogDensity(), stepSize);
            }
        }

        leapFrogEngine.updateMomentum(position, momentum,
                gradientProvider.getGradientLogDensity(), stepSize / 2);

        final double res = getScaledDotProduct(momentum, sigmaSquaredInverse) +
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
                            final double functionalStepSize,
                            final double[] sigmaSquaredInverse);

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
                                       double functionalStepSize, double[] sigmaSquaredInverse) {

                final int dim = momentum.length;
                for (int j = 0; j < dim; j++) {
                    position[j] += functionalStepSize * momentum[j] * sigmaSquaredInverse[j];
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
}
