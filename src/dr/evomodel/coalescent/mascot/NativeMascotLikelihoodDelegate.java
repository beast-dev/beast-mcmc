/*
 * NativeMascotLikelihoodDelegate.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import beagle.mascot.MascotJNIWrapper;
import dr.evomodel.coalescent.StructuredCoalescentSchedule;
import dr.evomodel.coalescent.StructuredCoalescentTipData;

/**
 * The MASCOT likelihood, gradients and ancestral state scores from beagle-lib's native component (see
 * {@link MascotJNIWrapper}). The native instance belongs to the {@link MascotLikelihood} for the whole run; this
 * delegate is a light wrapper around it with the current epoch boundaries, so it can be rebuilt as freely as
 * {@link GenericMascotLikelihoodDelegate}.
 * <p>
 * Migration matrices and inverse population sizes come from {@link MascotRuntime#updateEpochRates}, exactly as the
 * Java delegate computes them, so a strict native instance reproduces the Java gradients and ancestral state scores bit
 * for bit (logL can differ by an ulp of log(lambda) per coalescence, the platform's log). Stochastic mapping runs in
 * Java at {@code mappingMaxStep}.
 *
 * @author Marc A. Suchard
 */
public final class NativeMascotLikelihoodDelegate implements MascotLikelihoodDelegate {

    /** A native instance and its fixed dimensions; owned by one MascotLikelihood. */
    static final class Instance {
        final int id;
        final int stateCount;
        final int lineageCount;
        final int epochCount;
        final String details;
        private StructuredCoalescentTipData tipsSet;
        private boolean finalized;

        Instance(int stateCount, int lineageCount, int epochCount, boolean strict) {
            this.stateCount = stateCount;
            this.lineageCount = lineageCount;
            this.epochCount = epochCount;
            this.id = MascotJNIWrapper.createInstance(stateCount, lineageCount, epochCount,
                    strict ? MascotJNIWrapper.FP_STRICT : MascotJNIWrapper.FP_FMA);
            if (id < 0) {
                throw new IllegalStateException("native MASCOT instance creation failed with status " + id);
            }
            this.details = MascotJNIWrapper.getInstanceDetails(id);
        }

        void setTips(StructuredCoalescentTipData tipData) {
            if (tipData == tipsSet) {
                return;
            }
            double[] partials = new double[lineageCount * stateCount];
            for (int node = 0; node < lineageCount; node++) {
                if (tipData.hasTipPartials(node) || tipData.hasObservedState(node)) {
                    tipData.writeTipPartials(node, partials, node * stateCount);
                }
            }
            check(MascotJNIWrapper.setTipPartials(id, partials));
            tipsSet = tipData;
        }

        synchronized void release() {
            if (!finalized) {
                MascotJNIWrapper.finalizeInstance(id);
                finalized = true;
            }
        }

        @Override
        protected void finalize() throws Throwable {
            release();
            super.finalize();
        }
    }

    private final Instance instance;
    private final double[] boundaries;
    private final double maxStep;
    private final double mappingMaxStep;
    private final int parameterCount;
    private final MascotRuntime.EpochRates rates;
    private final double[] migration;
    private final double[] inversePopulation;
    private final double[] logLikelihood = new double[1];
    private GenericMascotLikelihoodDelegate mappingDelegate;

    NativeMascotLikelihoodDelegate(Instance instance, double[] boundaries, double maxStep, double mappingMaxStep) {
        if (boundaries.length - 1 != instance.epochCount) {
            throw new IllegalArgumentException("boundaries define " + (boundaries.length - 1) + " epochs, but the " +
                    "native instance holds " + instance.epochCount);
        }
        if (maxStep <= 0.0 || !Double.isFinite(maxStep)) {
            throw new IllegalArgumentException("maxStep must be positive and finite");
        }
        int k = instance.stateCount;
        this.instance = instance;
        this.boundaries = boundaries.clone();
        this.maxStep = maxStep;
        this.mappingMaxStep = mappingMaxStep;
        this.parameterCount = instance.epochCount * (k * (k - 1) + k);
        this.rates = new MascotRuntime.EpochRates(k);
        this.migration = new double[instance.epochCount * k * k];
        this.inversePopulation = new double[instance.epochCount * k];
    }

    @Override
    public double calculateLikelihood(MascotPreparedInput prepared, double[] theta, double[] branchRates,
                                      boolean checkProbabilities) {
        evaluate(prepared, theta, branchRates, checkProbabilities, MascotJNIWrapper.OUT_LOGL, null, null, null);
        return logLikelihood[0];
    }

    @Override
    public double calculateLikelihoodAndDerivatives(MascotPreparedInput prepared, double[] theta, double[] branchRates,
                                                    double[] gradientOut, double[] clockGradientOut,
                                                    double[] ancestralStateScoresOut, boolean checkProbabilities) {
        if (gradientOut != null && gradientOut.length != parameterCount) {
            throw new IllegalArgumentException("gradientOut dimension " + gradientOut.length +
                    " does not match expected dimension " + parameterCount);
        }
        if (clockGradientOut != null && branchRates == null) {
            throw new IllegalArgumentException("clockGradientOut was requested but branchRates is null");
        }
        int outputs = MascotJNIWrapper.OUT_LOGL |
                (gradientOut != null ? MascotJNIWrapper.OUT_THETA_GRADIENT : 0) |
                (clockGradientOut != null ? MascotJNIWrapper.OUT_CLOCK_GRADIENT : 0) |
                (ancestralStateScoresOut != null ? MascotJNIWrapper.OUT_ANCESTRAL_SCORES : 0);
        evaluate(prepared, theta, branchRates, checkProbabilities, outputs, gradientOut, clockGradientOut,
                ancestralStateScoresOut);
        return logLikelihood[0];
    }

    @Override
    public double prepareStochasticMap(MascotPreparedInput prepared, double[] theta, double[] branchRates,
                                       double[] nodeHeights, boolean checkProbabilities,
                                       StochasticMapSampler sampler) {
        if (mappingDelegate == null) {
            mappingDelegate = new GenericMascotLikelihoodDelegate(instance.stateCount, boundaries, mappingMaxStep);
        }
        return mappingDelegate.prepareStochasticMap(prepared, theta, branchRates, nodeHeights, checkProbabilities,
                sampler);
    }

    @Override
    public MascotStochasticMap drawStochasticMap() {
        if (mappingDelegate == null) {
            throw new IllegalStateException("prepareStochasticMap must succeed before drawStochasticMap");
        }
        return mappingDelegate.drawStochasticMap();
    }

    /**
     * A delegate with its own native instance, for tools and tests ({@link MascotLikelihood} shares one instance
     * across the delegates it rebuilds). Requires {@link MascotJNIWrapper#getUnavailableReason()} to be null.
     *
     * @param lineageCount the maximum lineage id + 1 of the inputs it will evaluate
     */
    public static NativeMascotLikelihoodDelegate create(int stateCount, int lineageCount, double[] boundaries,
                                                        double maxStep, boolean strict) {
        return new NativeMascotLikelihoodDelegate(new Instance(stateCount, lineageCount, boundaries.length - 1, strict),
                boundaries, maxStep, maxStep);
    }

    /** "kernel simdWidth fpMode isa", e.g. "A 0 fma neon". */
    public String getDetails() {
        return instance.details;
    }

    private void evaluate(MascotPreparedInput prepared, double[] theta, double[] branchRates,
                          boolean checkProbabilities, int outputs, double[] gradientOut, double[] clockGradientOut,
                          double[] ancestralStateScoresOut) {
        if (prepared.tipData.stateCount != instance.stateCount || prepared.maxLineageId + 1 != instance.lineageCount) {
            throw new IllegalArgumentException("prepared input does not match the native MASCOT instance");
        }
        if (theta == null || theta.length != parameterCount) {
            throw new IllegalArgumentException("theta dimension " + (theta == null ? -1 : theta.length) +
                    " does not match expected dimension " + parameterCount);
        }
        if (branchRates != null && branchRates.length < instance.lineageCount) {
            throw new IllegalArgumentException("branchRates dimension " + branchRates.length +
                    " does not cover the maximum lineage id " + prepared.maxLineageId);
        }
        int k = instance.stateCount;
        for (int epoch = 0; epoch < instance.epochCount; epoch++) {
            MascotRuntime.updateEpochRates(theta, epoch, k, rates);
            System.arraycopy(rates.migrationMatrix, 0, migration, epoch * k * k, k * k);
            System.arraycopy(rates.inversePopulation, 0, inversePopulation, epoch * k, k);
        }
        instance.setTips(prepared.tipData);
        StructuredCoalescentSchedule s = prepared.schedule;
        int status = MascotJNIWrapper.evaluate(instance.id, s.getIntervalCount(), s.initialTime,
                s.initialSampleNode, s.intervalLengths, s.eventTypes, s.sampleNodes, s.coalescentChild1,
                s.coalescentChild2, s.coalescentParent, instance.epochCount, boundaries, migration, inversePopulation,
                branchRates, MascotJNIWrapper.STEPS_FIXED, maxStep, 0.0, 0.0, 0, checkProbabilities, -1, -1, outputs,
                logLikelihood, gradientOut, clockGradientOut, ancestralStateScoresOut);
        check(status);
    }

    private static void check(int status) {
        switch (status) {
            case MascotJNIWrapper.OK:
                return;
            case MascotJNIWrapper.ERROR_COALESCENT_RATE:
                throw new NumericalException("invalid coalescent rate (native MASCOT)");
            case MascotJNIWrapper.ERROR_PROBABILITIES:
                throw new NumericalException("invalid lineage probabilities (native MASCOT)");
            default:
                throw new IllegalStateException("native MASCOT failed with status " + status);
        }
    }
}
