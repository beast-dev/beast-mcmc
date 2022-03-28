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
import dr.inference.hmc.PathGradient;
import dr.inference.hmc.ReversibleHMCProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.math.MathUtils;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.Transform;

/**
 * @author Max Tolkoff
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */

public class HamiltonianMonteCarloOperator extends AbstractAdaptableOperator
        implements GeneralOperator, PathDependent, ReversibleHMCProvider {

    final GradientWrtParameterProvider gradientProvider;
    protected double stepSize;
    LeapFrogEngine leapFrogEngine;
    protected final Parameter parameter;
    protected final MassPreconditioner preconditioning;
    protected final MassPreconditionScheduler preconditionScheduler;
    private final Options runtimeOptions;
    protected final double[] mask;
    protected final Transform transform;

//    public HamiltonianMonteCarloOperator(AdaptationMode mode, double weight,
//                                         GradientWrtParameterProvider gradientProvider,
//                                         Parameter parameter, Transform transform, Parameter maskParameter,
//                                         Options runtimeOptions,
//                                         MassPreconditioner.Type preconditioningType) {
//
//        super(mode, runtimeOptions.targetAcceptanceProbability);
//
//        setWeight(weight);
//
//        this.gradientProvider = gradientProvider;
//        this.runtimeOptions = runtimeOptions;
//        this.stepSize = runtimeOptions.initialStepSize;
//        this.preconditioning = preconditioningType.factory(gradientProvider, transform, runtimeOptions);
//        this.parameter = parameter;
//        this.mask = buildMask(maskParameter);
//        this.transform = transform;
//
//        this.leapFrogEngine = constructLeapFrogEngine(transform);
//    }
public HamiltonianMonteCarloOperator(AdaptationMode mode, double weight,
                                     GradientWrtParameterProvider gradientProvider,
                                     Parameter parameter, Transform transform, Parameter maskParameter,
                                     Options runtimeOptions,
                                     MassPreconditioner preconditioner) {
    this(mode, weight, gradientProvider, parameter, transform, maskParameter, runtimeOptions,
            preconditioner, MassPreconditionScheduler.Type.DEFAULT);
}

    public HamiltonianMonteCarloOperator(AdaptationMode mode, double weight,
                                         GradientWrtParameterProvider gradientProvider,
                                         Parameter parameter, Transform transform, Parameter maskParameter,
                                         Options runtimeOptions,
                                         MassPreconditioner preconditioner,
                                         MassPreconditionScheduler.Type preconditionSchedulerType) {

        super(mode, runtimeOptions.targetAcceptanceProbability);

        setWeight(weight);

        this.gradientProvider = gradientProvider;
        this.runtimeOptions = runtimeOptions;
        this.stepSize = runtimeOptions.initialStepSize;
        this.preconditioning = preconditioner;
        this.preconditionScheduler = preconditionSchedulerType.factory(runtimeOptions, (AdaptableMCMCOperator) this);
        this.parameter = parameter;
        this.mask = buildMask(maskParameter);
        this.transform = transform;

        this.leapFrogEngine = constructLeapFrogEngine(transform);
    }

    protected LeapFrogEngine constructLeapFrogEngine(Transform transform) {
        return (transform != null ?
                new LeapFrogEngine.WithTransform(parameter, transform,
                        getDefaultInstabilityHandler(), preconditioning, mask) :
                new LeapFrogEngine.Default(parameter,
                        getDefaultInstabilityHandler(), preconditioning, mask));
    }

    @Override
    public String getOperatorName() {
        return "VanillaHMC(" + parameter.getParameterName() + ")";
    }

    protected double[] buildMask(Parameter maskParameter) {

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

        if (preconditionScheduler.shouldUpdatePreconditioning()) {
            updatePreconditioning();
        }

        try {
            return leapFrog();
        } catch (NumericInstabilityException e) {
            return Double.NEGATIVE_INFINITY;
        } catch (ArithmeticException e) {
            if (REJECT_ARITHMETIC_EXCEPTION) {
                return Double.NEGATIVE_INFINITY;
            } else {
                throw e;
            }
        }
    }

    private void updatePreconditioning() {

        double[] lastGradient = leapFrogEngine.getLastGradient();
        double[] lastPosition = leapFrogEngine.getLastPosition();
        if (preconditionScheduler.shouldStoreSecant(lastGradient, lastPosition)) {
            preconditioning.storeSecant(new WrappedVector.Raw(lastGradient), new WrappedVector.Raw(lastPosition));
        }
        preconditioning.updateMass();
    }

    private static final boolean REJECT_ARITHMETIC_EXCEPTION = true;

    @Override
    public void setPathParameter(double beta) {
        if (gradientProvider instanceof PathGradient) {
            ((PathGradient) gradientProvider).setPathParameter(beta);
        }
    }

    private boolean shouldCheckStepSize() {
        return getCount() < 1 && getMode() == AdaptationMode.ADAPTATION_ON;
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
            } catch (Exception exception) {
                // Do nothing
            }

            if (!acceptableSize) {
                stepSize *= runtimeOptions.checkStepSizeReductionFactor;
            }

            ReadableVector.Utils.setParameter(initialPosition, parameter);  // Restore initial position
            ++iterations;
        }

        if (!acceptableSize && iterations < runtimeOptions.checkStepSizeMaxIterations) {
            throw new RuntimeException("Unable to find acceptable initial HMC step-size");
        }
    }

    boolean shouldCheckGradient() {
        return getCount() < runtimeOptions.gradientCheckCount;
    }

    void checkGradient(final Likelihood joint) {

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
                        "\tNumeric : " + new WrappedVector.Raw(numericGradientTransformed) + "\n" +
                        "\tParameter : " + new WrappedVector.Raw(parameter.getParameterValues()) + "\n" +
                        "\tTransformed Parameter : " + new WrappedVector.Raw(transformedParameter) + "\n";
                throw new RuntimeException(sb);
            }
        }

        ReadableVector.Utils.setParameter(restoredParameterValue, parameter);
    }

    static double[] mask(double[] vector, double[] mask) {

        assert (mask == null || mask.length == vector.length);

        if (mask != null) {
            for (int i = 0; i < vector.length; ++i) {
                if (mask[i] == 0.0) {
                    vector[i] = 0.0;
                }
            }
        }

        return vector;
    }

    static WrappedVector mask(WrappedVector vector, double[] mask) {

        assert (mask == null || mask.length == vector.getDim());

        if (mask != null) {
            for (int i = 0; i < vector.getDim(); ++i) {
                if (mask[i] == 0.0) {
                    vector.set(i, 0.0);
                }
            }
        }

        return vector;
    }

    private static final boolean DEBUG = false;

    public static class Options implements MassPreconditioningOptions {

        final double initialStepSize;
        final int nSteps;
        final double randomStepCountFraction;
        final int gradientCheckCount;
        final MassPreconditioningOptions preconditioningOptions;
        final double gradientCheckTolerance;
        final int checkStepSizeMaxIterations;
        final double checkStepSizeReductionFactor;
        final double targetAcceptanceProbability;
        final InstabilityHandler instabilityHandler;

        public Options(double initialStepSize, int nSteps, double randomStepCountFraction,
                       MassPreconditioningOptions preconditioningOptions,
                       int gradientCheckCount, double gradientCheckTolerance,
                       int checkStepSizeMaxIterations, double checkStepSizeReductionFactor,
                       double targetAcceptanceProbability, InstabilityHandler instabilityHandler) {
            this.initialStepSize = initialStepSize;
            this.nSteps = nSteps;
            this.randomStepCountFraction = randomStepCountFraction;
            this.gradientCheckCount = gradientCheckCount;
            this.gradientCheckTolerance = gradientCheckTolerance;
            this.checkStepSizeMaxIterations = checkStepSizeMaxIterations;
            this.checkStepSizeReductionFactor = checkStepSizeReductionFactor;
            this.targetAcceptanceProbability = targetAcceptanceProbability;
            this.instabilityHandler = instabilityHandler;
            this.preconditioningOptions = preconditioningOptions;
        }

        @Override
        public int preconditioningUpdateFrequency() {
            return preconditioningOptions.preconditioningUpdateFrequency();
        }

        @Override
        public int preconditioningDelay() {
            return preconditioningOptions.preconditioningDelay();
        }

        @Override
        public int preconditioningMaxUpdate() {
            return preconditioningOptions.preconditioningMaxUpdate();
        }

        @Override
        public int preconditioningMemory() {
            return preconditioningOptions.preconditioningMemory();
        }

        @Override
        public Parameter preconditioningEigenLowerBound() {
            throw new RuntimeException("Not yet implemented.");
        }

        @Override
        public Parameter preconditioningEigenUpperBound() {
            throw new RuntimeException("Not yet implemented.");
        }
    }

    public static class NumericInstabilityException extends Exception {
    }

    private int getNumberOfSteps() {
        int count = runtimeOptions.nSteps;
        if (runtimeOptions.randomStepCountFraction > 0.0) {
            double draw = count * (1.0 + runtimeOptions.randomStepCountFraction * (MathUtils.nextDouble() - 0.5));
            count = Math.max(1, (int) draw);
        }
        return count;
    }

    public double getKineticEnergy(ReadableVector momentum) {

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

        final WrappedVector momentum = mask(preconditioning.drawInitialMomentum(), mask);
        return leapFrogGivenMomentum(momentum);
    }

    protected double leapFrogGivenMomentum(WrappedVector momentum) throws NumericInstabilityException {
        leapFrogEngine.updateMask();
        final double[] position = leapFrogEngine.getInitialPosition();
        leapFrogEngine.projectMomentum(momentum.getBuffer(), position); //if momentum restricted to subspace

        final double prop = getKineticEnergy(momentum) +
                leapFrogEngine.getParameterLogJacobian();

        leapFrogEngine.updateMomentum(position, momentum.getBuffer(),
                mask(gradientProvider.getGradientLogDensity(), mask), stepSize / 2);


        int nStepsThisLeap = getNumberOfSteps();

        for (int i = 0; i < nStepsThisLeap; i++) { // Leap-frog

            try {
                leapFrogEngine.updatePosition(position, momentum, stepSize);
            } catch (ArithmeticException e) {
                throw new NumericInstabilityException();
            }

            if (i < (nStepsThisLeap - 1)) {

                try {
                    leapFrogEngine.updateMomentum(position, momentum.getBuffer(),
                            mask(gradientProvider.getGradientLogDensity(), mask), stepSize);
                } catch (ArithmeticException e) {
                    throw new NumericInstabilityException();
                }
            }
        }

        leapFrogEngine.updateMomentum(position, momentum.getBuffer(),
                mask(gradientProvider.getGradientLogDensity(), mask), stepSize / 2);

        final double res = getKineticEnergy(momentum) +
                leapFrogEngine.getParameterLogJacobian();

        return prop - res; //hasting ratio
    }

    @Override
    protected double getAdaptableParameterValue() {
        return Math.log(stepSize);
    }

    @Override
    public void setAdaptableParameterValue(double value) {
        if (DEBUG) {
            System.err.println("Setting adaptable parameter: " + getAdaptableParameter() + " -> " + value);
        }
        stepSize = Math.exp(value);
    }

    @Override
    public double getRawParameter() {
        return stepSize;
    }

    public enum InstabilityHandler {

        REJECT("reject") {
            @Override
            void checkValue(double x) throws NumericInstabilityException {
                if (Double.isNaN(x)) throw new NumericInstabilityException();
            }

            @Override
            void checkPosition(Transform transform, double[] unTransformedPosition) throws NumericInstabilityException {
                if (!transform.isInInteriorDomain(unTransformedPosition, 0, unTransformedPosition.length)) {
                    throw new NumericInstabilityException();
                }
            }

//            @Override
//            void checkEqual(double x, double y, double eps) throws NumericInstabilityException {
//                if (Math.abs(x - y) > eps) {
//                    throw new NumericInstabilityException();
//                }
//            }

            @Override
            boolean checkPositionTransform() {
                return true;
            }
        },

        DEBUG("debug") {
            @Override
            void checkValue(double x) throws NumericInstabilityException {
                if (Double.isNaN(x)) {
                    System.err.println("Numerical instability in HMC momentum; throwing exception");
                    throw new NumericInstabilityException();
                }
            }

            @Override
            void checkPosition(Transform transform, double[] unTransformedPosition) throws NumericInstabilityException {
                if (!transform.isInInteriorDomain(unTransformedPosition, 0, unTransformedPosition.length)) {
                    System.err.println("Numerical instability in HMC momentum; throwing exception");
                    throw new NumericInstabilityException();
                }
            }

//            @Override
//            void checkEqual(double x, double y, double eps) throws NumericInstabilityException {
//                if (Math.abs(x - y) > eps) {
//                    System.err.println("Numerical instability in HMC momentum; throwing exception");
//                    throw new NumericInstabilityException();
//                }
//            }

            @Override
            boolean checkPositionTransform() {
                return true;
            }
        },

        IGNORE("ignore") {
            @Override
            void checkValue(double x) {
                // Do nothing
            }

            @Override
            void checkPosition(Transform transform, double[] unTransformedPosition) throws NumericInstabilityException {
                // Do nothing
            }

//            @Override
//            void checkEqual(double x, double y, double eps) {
//                // Do nothing
//            }

            @Override
            boolean checkPositionTransform() {
                return false;
            }
        };

        private final String name;

        InstabilityHandler(String name) {
            this.name = name;
        }

        public static InstabilityHandler factory(String match) {
            for (InstabilityHandler type : InstabilityHandler.values()) {
                if (match.equalsIgnoreCase(type.name)) {
                    return type;
                }
            }
            return null;
        }

        abstract void checkValue(double x) throws NumericInstabilityException;

        //        abstract void checkEqual(double x, double y, double eps) throws NumericInstabilityException;
        abstract void checkPosition(Transform transform, double[] unTransformedPosition) throws NumericInstabilityException;

        abstract boolean checkPositionTransform();
    }

    protected InstabilityHandler getDefaultInstabilityHandler() {
        if (DEBUG) {
            return InstabilityHandler.DEBUG;
        } else {
            return runtimeOptions.instabilityHandler;
        }
    }

    @Override
    public String getAdaptableParameterName() {
        return "stepSize";
    }

    interface LeapFrogEngine {

        double[] getInitialPosition();

        double getParameterLogJacobian();

        void updateMomentum(final double[] position,
                            final double[] momentum,
                            final double[] gradient,
                            final double functionalStepSize) throws NumericInstabilityException;

        void updatePosition(final double[] position,
                            final WrappedVector momentum,
                            final double functionalStepSize) throws NumericInstabilityException;

        void setParameter(double[] position);

        double[] getLastGradient();

        double[] getLastPosition();

        void projectMomentum(double[] momentum, double[] position);

        void updateMask();

        class Default implements LeapFrogEngine {

            final protected Parameter parameter;
            final InstabilityHandler instabilityHandler;
            final private MassPreconditioner preconditioning;

            final double[] mask;

            double[] lastGradient;
            double[] lastPosition;

            Default(Parameter parameter, InstabilityHandler instabilityHandler,
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
            public void projectMomentum(double[] momentum, double[] position) {
                // do nothing
            }

            @Override
            public void updateMask() {
                // do nothing
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
            public void updatePosition(double[] position, WrappedVector momentum,
                                       double functionalStepSize) throws NumericInstabilityException {

                final int dim = momentum.getDim();
                for (int i = 0; i < dim; i++) {
                    position[i] += functionalStepSize * preconditioning.getVelocity(i, momentum);
                    instabilityHandler.checkValue(position[i]);
                }

                setParameter(position);
            }

            public void setParameter(double[] position) {
                ReadableVector.Utils.setParameter(position, parameter); // May not work with MaskedParameter?
            }
        }

        class WithTransform extends Default {

            final protected Transform transform;
            double[] unTransformedPosition;

            WithTransform(Parameter parameter, Transform transform,
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
            public void updatePosition(double[] position, WrappedVector momentum,
                                       double functionalStepSize) throws NumericInstabilityException {

                super.updatePosition(position, momentum, functionalStepSize);

                if (instabilityHandler.checkPositionTransform()) {
                    checkPosition(unTransformedPosition);
                }
            }

            @Override
            public void setParameter(double[] position) {
                unTransformedPosition = transform.inverse(position, 0, position.length);
                super.setParameter(unTransformedPosition);
            }

            private void checkPosition(double[] unTransformedPosition) throws NumericInstabilityException {
                instabilityHandler.checkPosition(transform, unTransformedPosition);
            }

//            private void checkPosition(double[] position) throws NumericInstabilityException {
//                double[] newPosition = transform.transform(transform.inverse(position, 0, position.length),
//                        0, position.length);
//                for (int i = 0; i < position.length; i++) {
//                    instabilityHandler.checkEqual(position[i], newPosition[i], EPS);
//                }
//            }
//
//            private double EPS = 10e-10;
        }
    }

    @Override
    public void reversiblePositionMomentumUpdate(WrappedVector position, WrappedVector momentum,
                                                 WrappedVector gradient, int direction, double time) {

        preconditionScheduler.forceUpdateCount();
        //providerUpdatePreconditioning();

        try {
            leapFrogEngine.updateMomentum(position.getBuffer(), momentum.getBuffer(),
                    mask(gradient.getBuffer(), mask), time * direction / 2);
            leapFrogEngine.updatePosition(position.getBuffer(), momentum, time * direction);
            updateGradient(gradient);
            leapFrogEngine.updateMomentum(position.getBuffer(), momentum.getBuffer(),
                    mask(gradient.getBuffer(), mask), time * direction / 2);
        } catch (NumericInstabilityException e) {
            handleInstability();
        }
    }

    @Override
    public void providerUpdatePreconditioning() {
        updatePreconditioning();
    }

    public void updateGradient(WrappedVector gradient) {
        double[] buffer = gradientProvider.getGradientLogDensity();
        for (int i = 0; i < buffer.length; i++) {
            gradient.set(i, buffer[i]);
        }
    }

    @Override
    public double[] getInitialPosition() {

        return leapFrogEngine.getInitialPosition();
    }

    @Override
    public double getParameterLogJacobian() {
        return leapFrogEngine.getParameterLogJacobian();
    }

    @Override
    public Transform getTransform() {
        return transform;
    }

    @Override
    public GradientWrtParameterProvider getGradientProvider() {
        return gradientProvider;
    }

    @Override
    public void setParameter(double[] position) {
        leapFrogEngine.setParameter(position);
    }

    @Override
    public WrappedVector drawMomentum() {
        return mask(preconditioning.drawInitialMomentum(), mask);
    }

    @Override
    public double getJointProbability(WrappedVector momentum) {
        return gradientProvider.getLikelihood().getLogLikelihood() - getKineticEnergy(momentum) - getParameterLogJacobian();
    }

    @Override
    public double getLogLikelihood() {
        return gradientProvider.getLikelihood().getLogLikelihood();
    }

    @Override
    public double getStepSize() {
        return stepSize;
    }

    public int getNumGradientEvent(){
        return 0;
    }

    @Override
    public int getNumBoundaryEvent() {
        return 0;
    }

    @Override
    public double[] getMask() {
        return mask;
    }

    protected void handleInstability() {
        throw new RuntimeException("Numerical instability; need to handle"); // TODO
    }
}
