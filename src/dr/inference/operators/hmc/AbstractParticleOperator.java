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

import dr.evomodel.operators.NativeZigZagOptions;
import dr.evomodel.operators.NativeZigZagWrapper;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.PrecisionColumnProvider;
import dr.inference.hmc.PrecisionMatrixVectorProductProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.BenchmarkTimer;
import dr.xml.Reportable;

import java.util.Arrays;

import static dr.math.matrixAlgebra.ReadableVector.Utils.setParameter;

/**
 * @author Aki Nishimura
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */

public abstract class AbstractParticleOperator extends SimpleMCMCOperator implements GibbsOperator, Reportable {

    private static final boolean CHECK_MATRIX_ILL_CONDITIONED = false;

    AbstractParticleOperator(GradientWrtParameterProvider gradientProvider,
                             PrecisionMatrixVectorProductProvider multiplicationProvider,
                             PrecisionColumnProvider columnProvider,
                             double weight, Options runtimeOptions, NativeCodeOptions nativeOptions,
                             boolean refreshVelocity, Parameter mask, Parameter categoryClass,
                             MassPreconditioner massPreconditioner,
                             MassPreconditionScheduler.Type preconditionSchedulerType) {

        this.gradientProvider = gradientProvider;
        this.productProvider = multiplicationProvider;
        this.columnProvider = columnProvider;
        this.parameter = gradientProvider.getParameter();
        this.mask = mask;
        this.maskVector = mask != null ? mask.getParameterValues() : null;
        this.parameterSign = setParameterSign(gradientProvider);

        this.runtimeOptions = runtimeOptions;
        this.nativeCodeOptions = nativeOptions;
        this.refreshVelocity = refreshVelocity;
        this.preconditioning = setupPreconditioning();
        this.meanVector = getMeanVector(gradientProvider);

        this.massPreconditioning = massPreconditioner;
        this.preconditionScheduler = preconditionSchedulerType.factory(runtimeOptions, this);

        setWeight(weight);
        this.observedDataMask = getObservedDataMask();
        this.categoryClasses = getCategoryClasses(categoryClass);
        checkParameterBounds(parameter);

        long flags = 128;
        long nativeSeed = MathUtils.nextLong();
        int nThreads = 4;

        if (nativeOptions.testNativeFindNextBounce || nativeOptions.useNativeFindNextBounce || nativeOptions.useNativeUpdateDynamics) {

            NativeZigZagOptions options = new NativeZigZagOptions(flags, nativeSeed, nThreads);

            nativeZigZag = new NativeZigZagWrapper(parameter.getDimension(), options,
                    maskVector, getObservedDataMask(), parameterSign);
        }
    }

    private double[] setParameterSign(GradientWrtParameterProvider gradientProvider) {

        double[] startingValue = gradientProvider.getParameter().getParameterValues();
        double[] sign = new double[startingValue.length];

        for (int i = 0; i < startingValue.length; i++) {

            if (startingValue[i] == 0 && (mask == null || mask.getParameterValue(i) == 1)) {
                throw new RuntimeException("Must start from either positive or negative value!");
            }
            sign[i] = startingValue[i] > 0 ? 1 : -1;
        }
        return sign;
    }

    private double[] getObservedDataMask() {
        int dim = parameter.getDimension();
        double[] observed = new double[dim];
        assert (dim == parameter.getBounds().getBoundsDimension());

        for (int i = 0; i < dim; ++i) {
            observed[i] = (parameter.getBounds().getUpperLimit(i) == Double.POSITIVE_INFINITY &&
                    parameter.getBounds().getLowerLimit(i) == Double.NEGATIVE_INFINITY) ? 0.0 : 1.0;
        }
        return observed;
    }

    private int[] getCategoryClasses(Parameter categoryVector) {
        int dim = parameter.getDimension();
        int[] categoryClasses = new int[dim];

        if (categoryVector != null) {
            int[] category = new int[categoryVector.getDimension()];
            for (int i = 0; i < category.length; i++) {
                category[i] = (int) categoryVector.getParameterValues()[i];
            }

            int L = categoryVector.getDimension();
            int n = dim / L;
            for (int i = 0; i < n; i++) {
                System.arraycopy(category, 0, categoryClasses, i * L, L);
            }
        }
        return categoryClasses;
    }

    @Override
    public double doOperation() {

        WrappedVector position = getInitialPosition();

        WrappedVector momentum = drawInitialMomentum();

        if (preconditionScheduler.shouldUpdatePreconditioning()){
            updatePreconditioning(position);
        }

        double hastingsRatio = integrateTrajectory(position, momentum);

        setParameter(position, parameter);

        if (CHECK_MATRIX_ILL_CONDITIONED & getCount() % 100 == 0) {
            productProvider.getTimeScaleEigen();
        }

        return hastingsRatio;
    }

    abstract double integrateTrajectory(WrappedVector position, WrappedVector momentum);

    WrappedVector drawInitialMomentum() {
        return new WrappedVector.Raw(null, 0, 0);
    }

    double drawTotalTravelTime() {
        double randomFraction = 1.0 + runtimeOptions.randomTimeWidth * (MathUtils.nextDouble() - 0.5);
        return preconditioning.totalTravelTime * randomFraction;
    }

    static void updateGradient(WrappedVector gradient, double time, WrappedVector action) {

        final double[] g = gradient.getBuffer();
        final double[] a = action.getBuffer();

        for (int i = 0, len = g.length; i < len; ++i) {
            g[i] -= time * a[i];
        }
    }

    static void updatePosition(WrappedVector position, WrappedVector velocity, double time) {

        final double[] p = position.getBuffer();
        final double[] v = velocity.getBuffer();

        for (int i = 0, len = p.length; i < len; ++i) {
            p[i] += time * v[i];
        }
    }

    static void updatePosition(double[] p, double[] v, double time) {
        for (int i = 0, len = p.length; i < len; ++i) {
            p[i] += time * v[i];
        }
    }


    static void updateMomentum(double[] a, double[] g, double[] m, double time) {

        final double halfTimeSquared = time * time / 2;

        for (int i = 0, len = m.length; i < len; ++i) {
            m[i] = m[i] + time * g[i] - halfTimeSquared * a[i];
        }
    }


    WrappedVector getInitialGradient() {

        double[] gradient = gradientProvider.getGradientLogDensity();

        if (mask != null) {
            applyMask(gradient);
        }

        return new WrappedVector.Raw(gradient);
    }

    @SuppressWarnings("unused")
    void applyMask(WrappedVector vector) {
        applyMask(vector.getBuffer());
    }

    void applyMask(double[] vector) {

        if (TIMING) {
            timer.startTimer("applyMask");
        }
        assert (vector.length == mask.getDimension());

        for (int i = 0, len = vector.length; i < len; ++i) {
            vector[i] *= maskVector[i];
        }

        if (TIMING) {
            timer.stopTimer("applyMask");
        }
    }

    WrappedVector getPrecisionProduct(ReadableVector velocity) {

        WrappedVector velocityPlusMean = new WrappedVector.Raw(new double[velocity.getDim()]);
        for (int i = 0; i < velocityPlusMean.getDim(); i++) {
            velocityPlusMean.set(i, velocity.get(i) + meanVector[i]);
        }
        setParameter(velocityPlusMean, parameter);
        double[] product = productProvider.getProduct(parameter);

        if (mask != null) {
            applyMask(product);
        }

        return new WrappedVector.Raw(product);
    }

    WrappedVector getPrecisionColumn(int index) {

        if (TIMING) {
            timer.startTimer("getColumn");
        }

        double[] precisionColumn = columnProvider.getColumn(index);

        if (TIMING) {
            timer.stopTimer("getColumn");
        }

        if (mask != null) {
            applyMask(precisionColumn);
        }

        return new WrappedVector.Raw(precisionColumn);
    }

    void updateAction(WrappedVector action, ReadableVector velocity, int eventIndex) {

        WrappedVector column = getPrecisionColumn(eventIndex);

        if (TIMING) {
            timer.startTimer("updateAction");
        }

        final double[] a = action.getBuffer();
        final double[] c = column.getBuffer();

        final double twoV = 2 * velocity.get(eventIndex);

        for (int i = 0, len = a.length; i < len; ++i) {
            a[i] += twoV * c[i];
        }

        if (TIMING) {
            timer.stopTimer("updateAction");
        }

        if (mask != null) {
            applyMask(a);
        }
    }

    boolean headingTowardsBinaryBoundary(double velocity, int positionIndex) {
        return observedDataMask[positionIndex] * parameterSign[positionIndex] * velocity < 0.0;
    }

    private WrappedVector getInitialPosition() {
        return new WrappedVector.Raw(parameter.getParameterValues());
    }

    private void checkParameterBounds(Parameter parameter) {

        for (int i = 0, len = parameter.getDimension(); i < len; ++i) {
            double value = parameter.getParameterValue(i);
            if (value < parameter.getBounds().getLowerLimit(i) ||
                    value > parameter.getBounds().getUpperLimit(i)) {
                throw new IllegalArgumentException("Parameter '" + parameter.getId() + "' is out-of-bounds");
            }
        }
    }

    private Preconditioning setupPreconditioning() {

        double[] mass = new double[parameter.getDimension()];
        Arrays.fill(mass, 1.0);

        // TODO Should use:
        productProvider.getMassVector();
        double time = productProvider.getTimeScale();

        return new Preconditioning(
                new WrappedVector.Raw(mass),
                time
        );
    }

    void updatePreconditioning(WrappedVector position) {
        massPreconditioning.updateVariance(position);
        massPreconditioning.updateMass();
        preconditioning.mass = massPreconditioning.getMass();
    }

    void initializeNumEvent() {
        numEvents = 0;
        numBoundaryEvents = 0;
        numGradientEvents = 0;
    }

    void recordOneMoreEvent() {
        numEvents++;
    }

    void recordEvents(Type eventType) {
        numEvents++;
        if (eventType == Type.BINARY_BOUNDARY || eventType == Type.CATE_BOUNDARY) {
            numBoundaryEvents++;
        } else if (eventType == Type.GRADIENT) {
            numGradientEvents++;
        }
    }

    void storeVelocity(WrappedVector velocity) {
        storedVelocity = velocity;
    }

    double[] getMeanVector(GradientWrtParameterProvider gradientProvider) {

        double[] mean = new double[parameter.getDimension()];

        if (gradientProvider.getLikelihood() instanceof TreeDataLikelihood) {

            TreeDataLikelihood likelihood = (TreeDataLikelihood) gradientProvider.getLikelihood();
            ContinuousDataLikelihoodDelegate likelihoodDelegate =
                    (ContinuousDataLikelihoodDelegate) likelihood.getDataLikelihoodDelegate();
            double[] rootMean = likelihoodDelegate.getRootPrior().getMean();
            int dimTrait = likelihoodDelegate.getTraitDim();
            int taxonCount = parameter.getDimension() / dimTrait;

            int index = 0;
            for (int taxon = 0; taxon < taxonCount; ++taxon) {
                for (int trait = 0; trait < dimTrait; ++trait) {
                    mean[index + trait] = rootMean[trait];
                }
                index += dimTrait;
            }
        }
        return mean;
    }

    public static class Options implements MassPreconditioningOptions {

        final double randomTimeWidth;
        final int preconditioningUpdateFrequency;
        final int preconditioningMaxUpdate;
        final int preconditioningDelay;
        final int updateSampleCovFrequency;
        final int updateSampleCovDelay;

        public Options(double randomTimeWidth, int preconditioningUpdateFrequency, int preconditioningMaxUpdate,
                       int preconditioningDelay, int updateSampleCovFrequency, int updateSampleCovDelay) {
            this.randomTimeWidth = randomTimeWidth;
            this.preconditioningUpdateFrequency = preconditioningUpdateFrequency;
            this.preconditioningMaxUpdate = preconditioningMaxUpdate;
            this.preconditioningDelay = preconditioningDelay;
            this.updateSampleCovFrequency = updateSampleCovFrequency;
            this.updateSampleCovDelay = updateSampleCovDelay;
        }

        @Override
        public int preconditioningUpdateFrequency() {
            return preconditioningUpdateFrequency;
        }

        @Override
        public int preconditioningDelay() {
            return preconditioningDelay;
        }

        @Override
        public int preconditioningMaxUpdate() {
            return preconditioningMaxUpdate;
        }

        @Override
        public int preconditioningMemory() {
            return 0;
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

    public static class NativeCodeOptions {
        final boolean testNativeFindNextBounce;
        final boolean useNativeFindNextBounce;
        final boolean useNativeUpdateDynamics;

        public NativeCodeOptions(boolean testNativeFindNextBounce, boolean useNativeFindNextBounce,
                                 boolean useNativeUpdateDynamics) {
            this.testNativeFindNextBounce = testNativeFindNextBounce;
            this.useNativeFindNextBounce = useNativeFindNextBounce;
            this.useNativeUpdateDynamics = useNativeUpdateDynamics;
        }
    }

    protected class Preconditioning {

        WrappedVector mass;
        double totalTravelTime;

        private Preconditioning(WrappedVector mass, double totalTravelTime) {
            this.mass = mass;
            this.totalTravelTime = totalTravelTime;
        }
    }

    class BounceState {
        final Type type;
        final int index;
        final double remainingTime;

        BounceState(Type type, int index, double remainingTime) {
            this.type = type;
            this.index = index;
            this.remainingTime = remainingTime;
        }

        BounceState(double remainingTime) {
            this.type = Type.NONE;
            this.index = -1;
            this.remainingTime = remainingTime;
        }

        boolean isTimeRemaining() {
            return remainingTime > 0.0;
        }

        public String toString() {
            return "remainingTime : " + remainingTime +
                    " lastBounceType: " + type + " in dim: " + index;
        }
    }

    enum Type {
        NONE,
        BINARY_BOUNDARY,
        CATE_BOUNDARY,
        GRADIENT,
        REFRESHMENT;

        public static Type castFromInt(int i) {
            if (i == 0) {
                return NONE;
            } else if (i == 1) {
                return BINARY_BOUNDARY;
            } else if (i == 2) {
                return GRADIENT;
            } else if (i > 2) {
                return CATE_BOUNDARY;
            } else {
                throw new RuntimeException("Unknown type");
            }
        }
    }

    @Override
    public String getReport() {
        return TIMING ? timer.toString() : "";
    }

    protected final GradientWrtParameterProvider gradientProvider;
    private final PrecisionMatrixVectorProductProvider productProvider;
    private final PrecisionColumnProvider columnProvider;
    protected final Parameter parameter;
    protected final Options runtimeOptions;
    protected boolean refreshVelocity;
    protected final NativeCodeOptions nativeCodeOptions;
    final Parameter mask;
    final double[] parameterSign;
    protected final double[] maskVector;
    int numEvents;
    int numBoundaryEvents;
    int numGradientEvents;
    protected WrappedVector storedVelocity;
    Preconditioning preconditioning;
    protected final MassPreconditioner massPreconditioning;
    protected final MassPreconditionScheduler preconditionScheduler;
    final protected double[] observedDataMask;
    private final double[] meanVector;

    protected final int[] categoryClasses;
    final static boolean TIMING = true;
    BenchmarkTimer timer = new BenchmarkTimer();

//  final static boolean TEST_CRITICAL_REGION = false;
    NativeZigZagWrapper nativeZigZag;
}
