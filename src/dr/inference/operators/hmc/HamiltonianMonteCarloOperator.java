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

import dr.inference.hmc.PathGradient;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.GeneralOperator;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;
import dr.inference.operators.PathDependent;
import dr.math.MathUtils;
import dr.util.Transform;

/**
 * @author Max Tolkoff
 * @author Marc A. Suchard
 */

public class HamiltonianMonteCarloOperator extends AbstractCoercableOperator
        implements GeneralOperator, PathDependent {

    final GradientWrtParameterProvider gradientProvider;
    protected double stepSize;
    final LeapFrogEngine leapFrogEngine;
    private final Parameter parameter;
    final MassPreconditioner preconditioning;
    private final Options runtimeOptions;
    protected final double[] mask;
    protected final Transform transform;

    HamiltonianMonteCarloOperator(CoercionMode mode, double weight, GradientWrtParameterProvider gradientProvider,
                                  Parameter parameter, Transform transform, Parameter mask,
                                  double stepSize, int nSteps,
                                  double randomStepCountFraction,
                                  double gradientCheckTolerance) {
        this(mode, weight, gradientProvider,
                parameter, transform, mask,
                new Options(stepSize, nSteps, randomStepCountFraction,
                        0, 0, 0, gradientCheckTolerance,
                        10, 0.1),
                MassPreconditioner.Type.NONE
        );
    }

    public HamiltonianMonteCarloOperator(CoercionMode mode, double weight,
                                         GradientWrtParameterProvider gradientProvider,
                                         Parameter parameter, Transform transform, Parameter maskParameter,
                                         Options runtimeOptions,
                                         MassPreconditioner.Type preconditioningType) {

        super(mode);
        setWeight(weight);
        setTargetAcceptanceProbability(0.8); // Stan default

        this.gradientProvider = gradientProvider;
        this.runtimeOptions = runtimeOptions;
        this.stepSize = runtimeOptions.initialStepSize;
        this.preconditioning = preconditioningType.factory(gradientProvider, transform);
        this.parameter = parameter;
        this.mask = buildMask(maskParameter);
        this.transform = transform;

        this.leapFrogEngine = (transform != null ?
                new LeapFrogEngine.WithTransform(parameter, transform,
                        getDefaultInstabilityHandler(), preconditioning, mask) :
                new LeapFrogEngine.Default(parameter,
                        getDefaultInstabilityHandler(), preconditioning, mask));
    }

    @Override
    public String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public String getOperatorName() {
        return "Vanilla HMC operator";
    }

    private boolean shouldUpdatePreconditioning() {
        return ((runtimeOptions.preconditioningUpdateFrequency > 0)
                && (((getCount() % runtimeOptions.preconditioningUpdateFrequency == 0)
                && (getCount() > runtimeOptions.preconditioningDelay))
                || (getCount() == runtimeOptions.preconditioningDelay)));
    }

    private static double[] buildMask(Parameter maskParameter) {

        if (maskParameter == null) return null;

        double[] mask = new double[maskParameter.getDimension()];

        for (int i = 0; i < mask.length; ++i) {
            mask[i] = (maskParameter.getParameterValue(i) == 0.0) ? 0.0 : 1.0;
        }

        return mask;
    }

    @Override
    public double doOperation() {
        throw new RuntimeException("Should not be executed");
    }

    @Override
    public double doOperation(Likelihood joint) {

        if (shouldCheckStepSize()) {
            checkStepSize();
        }

        if (shouldCheckGradient()) {
            checkGradient(joint);

        }

        if (shouldUpdatePreconditioning()) {
            preconditioning.updateMass();
        }

        try {
            return leapFrog();
        } catch (NumericInstabilityException e) {
            return Double.NEGATIVE_INFINITY;
        }
    }

    @Override
    public void setPathParameter(double beta) {
        if (gradientProvider instanceof PathGradient) {
            ((PathGradient) gradientProvider).setPathParameter(beta);
        }
    }

    private boolean shouldCheckStepSize() {
        return getCount() < 1 && getMode() == CoercionMode.COERCION_ON;
    }

    private void checkStepSize() {

        double[] initialPosition = parameter.getParameterValues();

        int iterations = 0;
        boolean acceptableSize = false;

        while (!acceptableSize && iterations < runtimeOptions.checkStepSizeMaxIterations) {

            try {
                leapFrog();
                double logLikelihood = gradientProvider.getLikelihood().getLogLikelihood();

                if (!Double.isNaN(logLikelihood) && !Double.isInfinite(logLikelihood)) {
                    acceptableSize = true;
                }
            } catch (AssertionError error) {
                // Do nothing
            } catch (Exception exception) {
                // Do nothing
            }

            if (!acceptableSize) {
                stepSize *= runtimeOptions.checkStepSizeMaxIterations;
            }

            ReadableVector.Utils.setParameter(initialPosition, parameter);  // Restore initial position
            ++iterations;
        }

        if (!acceptableSize) {
            throw new RuntimeException("Unable to find acceptable initial HMC step-size");
        }
    }

    private boolean shouldCheckGradient() {
        return getCount() < runtimeOptions.gradientCheckCount;
    }

    private void checkGradient(final Likelihood joint) {

        if (parameter.getDimension() != gradientProvider.getDimension()) {
            throw new RuntimeException("Unequal dimensions");
        }

        MultivariateFunction numeric = new MultivariateFunction() {

            @Override
            public double evaluate(double[] argument) {

                if (transform == null) {

                    ReadableVector.Utils.setParameter(argument, parameter);
                    return joint.getLogLikelihood();
                } else {

                    double[] untransformedValue = transform.inverse(argument, 0, argument.length);
                    ReadableVector.Utils.setParameter(untransformedValue, parameter);
                    return joint.getLogLikelihood() - transform.getLogJacobian(untransformedValue, 0, untransformedValue.length);
                }
            }

            @Override
            public int getNumArguments() {
                return parameter.getDimension();
            }

            @Override
            public double getLowerBound(int n) {
                return parameter.getBounds().getLowerLimit(n);
            }

            @Override
            public double getUpperBound(int n) {
                return parameter.getBounds().getUpperLimit(n);
            }
        };

        double[] analyticalGradientOriginal = gradientProvider.getGradientLogDensity();
        double[] restoredParameterValue = parameter.getParameterValues();

        if (transform == null) {

            double[] numericGradientOriginal = NumericalDerivative.gradient(numeric, parameter.getParameterValues());

            if (!MathUtils.isClose(analyticalGradientOriginal, numericGradientOriginal, runtimeOptions.gradientCheckTolerance)) {

                String sb = "Gradients do not match:\n" +
                        "\tAnalytic: " + new WrappedVector.Raw(analyticalGradientOriginal) + "\n" +
                        "\tNumeric : " + new WrappedVector.Raw(numericGradientOriginal) + "\n";
                throw new RuntimeException(sb);
            }

        } else {

            double[] transformedParameter = transform.transform(parameter.getParameterValues(), 0,
                    parameter.getParameterValues().length);
            double[] numericGradientTransformed = NumericalDerivative.gradient(numeric, transformedParameter);

            double[] analyticalGradientTransformed = transform.updateGradientLogDensity(analyticalGradientOriginal,
                    parameter.getParameterValues(), 0, parameter.getParameterValues().length);

            if (!MathUtils.isClose(analyticalGradientTransformed, numericGradientTransformed, runtimeOptions.gradientCheckTolerance)) {
                String sb = "Transformed Gradients do not match:\n" +
                        "\tAnalytic: " + new WrappedVector.Raw(analyticalGradientTransformed) + "\n" +
                        "\tNumeric : " + new WrappedVector.Raw(numericGradientTransformed) + "\n";
                throw new RuntimeException(sb);
            }
        }

        ReadableVector.Utils.setParameter(restoredParameterValue, parameter);
    }


    static double[] mask(double[] vector, double[] mask) {

        assert (mask == null || mask.length == vector.length);

        if (mask != null) {
            for (int i = 0; i < vector.length; ++i) {
                vector[i] *= mask[i];
            }
        }

        return vector;
    }

    static WrappedVector mask(WrappedVector vector, double[] mask) {

        assert (mask == null || mask.length == vector.getDim());

        if (mask != null) {
            for (int i = 0; i < vector.getDim(); ++i) {
                vector.set(i, vector.get(i) * mask[i]);
            }
        }

        return vector;
    }

    private static final boolean DEBUG = false;

    public static class Options {

        final double initialStepSize;
        final int nSteps;
        final double randomStepCountFraction;
        final int preconditioningUpdateFrequency;
        final int preconditioningDelay;
        final int gradientCheckCount;
        final double gradientCheckTolerance;
        final int checkStepSizeMaxIterations;
        final double checkStepSizeReductionFactor;

        public Options(double initialStepSize, int nSteps, double randomStepCountFraction, int preconditioningUpdateFrequency,
                       int preconditioningDelay, int gradientCheckCount, double gradientCheckTolerance,
                       int checkStepSizeMaxIterations, double checkStepSizeReductionFactor) {
            this.initialStepSize = initialStepSize;
            this.nSteps = nSteps;
            this.randomStepCountFraction = randomStepCountFraction;
            this.preconditioningUpdateFrequency = preconditioningUpdateFrequency;
            this.preconditioningDelay = preconditioningDelay;
            this.gradientCheckCount = gradientCheckCount;
            this.gradientCheckTolerance = gradientCheckTolerance;
            this.checkStepSizeMaxIterations = checkStepSizeMaxIterations;
            this.checkStepSizeReductionFactor = checkStepSizeReductionFactor;
        }
    }

    static class NumericInstabilityException extends Exception { }

    private int getNumberOfSteps() {
        int count = runtimeOptions.nSteps;
        if (runtimeOptions.randomStepCountFraction > 0.0) {
            double draw = count * (1.0 + runtimeOptions.randomStepCountFraction * (MathUtils.nextDouble() - 0.5));
            count = Math.max(1, (int) draw);
        }
        return count;
    }

    double getKineticEnergy(ReadableVector momentum) {

        final int dim = momentum.getDim();

        double energy = 0.0;
        for (int i = 0; i < dim; i++) {
            energy += momentum.get(i) * preconditioning.getVelocity(i, momentum);
        }
        return energy / 2.0;
    }

    private double leapFrog() throws NumericInstabilityException {

        if (DEBUG) {
            System.err.println("HMC step size: " + stepSize);
        }

        final double[] position = leapFrogEngine.getInitialPosition();
        final WrappedVector momentum = mask(preconditioning.drawInitialMomentum(), mask);

        final double prop = getKineticEnergy(momentum) +
                leapFrogEngine.getParameterLogJacobian();

        leapFrogEngine.updateMomentum(position, momentum.getBuffer(),
                mask(gradientProvider.getGradientLogDensity(), mask), stepSize / 2);


        int nStepsThisLeap = getNumberOfSteps();

        for (int i = 0; i < nStepsThisLeap; i++) { // Leap-frog

            leapFrogEngine.updatePosition(position, momentum, stepSize);

            if (i < (nStepsThisLeap - 1)) {
                leapFrogEngine.updateMomentum(position, momentum.getBuffer(),
                        mask(gradientProvider.getGradientLogDensity(), mask), stepSize);
            }
        }

        leapFrogEngine.updateMomentum(position, momentum.getBuffer(),
                mask(gradientProvider.getGradientLogDensity(), mask), stepSize / 2);

        final double res = getKineticEnergy(momentum) +
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

    @Override
    public void accept(double deviation) {

        super.accept(deviation);
        preconditioning.storeSecant(
                new WrappedVector.Raw(leapFrogEngine.getLastGradient()),
                new WrappedVector.Raw(leapFrogEngine.getLastPosition())
        );
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
                            final ReadableVector momentum,
                            final double functionalStepSize);

        void setParameter(double[] position);

        double[] getLastGradient();

        double[] getLastPosition();

        class Default implements LeapFrogEngine {

            final protected Parameter parameter;
            final private InstabilityHandler instabilityHandler;
            final private MassPreconditioner preconditioning;

            final double[] mask;

            double[] lastGradient;
            double[] lastPosition;

            protected Default(Parameter parameter, InstabilityHandler instabilityHandler,
                              MassPreconditioner preconditioning,
                              double[] mask) {
                this.parameter = parameter;
                this.instabilityHandler = instabilityHandler;
                this.preconditioning = preconditioning;
                this.mask = mask;
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
            public double[] getLastGradient() {
                return lastGradient;
            }

            @Override
            public double[] getLastPosition() {
                return lastPosition;
            }

            @Override
            public void updateMomentum(double[] position, double[] momentum, double[] gradient,
                                       double functionalStepSize) throws NumericInstabilityException {

                final int dim = momentum.length;
                for (int i = 0; i < dim; ++i) {
                    momentum[i] += functionalStepSize * gradient[i];
                    instabilityHandler.checkValue(momentum[i]);
                }

                lastGradient = gradient;
                lastPosition = position;
            }

            @Override
            public void updatePosition(double[] position, ReadableVector momentum,
                                       double functionalStepSize) {

                final int dim = momentum.getDim();
                for (int i = 0; i < dim; i++) {
                    position[i] += functionalStepSize * preconditioning.getVelocity(i, momentum);
                }
                setParameter(position);
            }

            public void setParameter(double[] position) {
                ReadableVector.Utils.setParameter(position, parameter); // May not work with MaskedParameter?
            }
        }

        class WithTransform extends Default {

            final private Transform transform;
            double[] unTransformedPosition;

            private WithTransform(Parameter parameter, Transform transform,
                                  InstabilityHandler instabilityHandler,
                                  MassPreconditioner preconditioning,
                                  double[] mask) {
                super(parameter, instabilityHandler, preconditioning, mask);
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
                                       double functionalStepSize) throws NumericInstabilityException {

                gradient = transform.updateGradientLogDensity(gradient, unTransformedPosition,
                        0, unTransformedPosition.length);
                mask(gradient, mask);
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
