/*
 * GenericMascotLikelihoodDelegate.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 */

package dr.evomodel.coalescent.mascot;

import java.util.Arrays;

/**
 * Allocation-light Java MASCOT likelihood delegate. It coordinates forward
 * evaluation, epoch-rate caching, and optional adjoint replay.
 */
public final class GenericMascotLikelihoodDelegate implements MascotLikelihoodDelegate {

    private final int stateCount;
    private final double[] boundaries;
    private final int migrationParametersPerEpoch;
    private final int parametersPerEpoch;
    private final int parameterCount;
    private final int epochCount;

    private final MascotRuntime.Workspace workspace = new MascotRuntime.Workspace();
    private final MascotRuntime.EpochRates[] epochRates;
    private final MascotRuntime.ActiveState activeState = new MascotRuntime.ActiveState();
    private final MascotForwardEvaluator forwardEvaluator;
    private final MascotAdjointEvaluator adjointEvaluator;

    public GenericMascotLikelihoodDelegate(int stateCount, double[] boundaries, double maxStep) {
        if (stateCount < 2) {
            throw new IllegalArgumentException("stateCount must be at least 2");
        }
        if (boundaries == null || boundaries.length < 2) {
            throw new IllegalArgumentException("boundaries must include at least 0.0 and infinity");
        }
        if (Math.abs(boundaries[0]) > MascotRuntime.TIME_TOLERANCE) {
            throw new IllegalArgumentException("the first boundary must be 0.0");
        }
        for (int i = 1; i < boundaries.length; i++) {
            if (!(boundaries[i] > boundaries[i - 1])) {
                throw new IllegalArgumentException("boundaries must be strictly increasing");
            }
        }
        if (!Double.isInfinite(boundaries[boundaries.length - 1])) {
            throw new IllegalArgumentException("the final boundary must be positive infinity");
        }
        if (maxStep <= 0.0 || !Double.isFinite(maxStep)) {
            throw new IllegalArgumentException("maxStep must be positive and finite");
        }

        this.stateCount = stateCount;
        this.boundaries = boundaries.clone();
        this.migrationParametersPerEpoch = stateCount * (stateCount - 1);
        this.parametersPerEpoch = migrationParametersPerEpoch + stateCount;
        this.epochCount = this.boundaries.length - 1;
        this.parameterCount = epochCount * parametersPerEpoch;

        this.epochRates = new MascotRuntime.EpochRates[epochCount];
        for (int epoch = 0; epoch < epochCount; epoch++) {
            epochRates[epoch] = new MascotRuntime.EpochRates(stateCount);
        }

        this.forwardEvaluator = new MascotForwardEvaluator(stateCount, maxStep, this.boundaries, epochCount);
        this.adjointEvaluator = new MascotAdjointEvaluator(stateCount, parametersPerEpoch,
                migrationParametersPerEpoch, parameterCount);
    }

    public int getStateCount() {
        return stateCount;
    }

    public int getEpochCount() {
        return epochCount;
    }

    public int getParameterCount() {
        return parameterCount;
    }

    public double logLikelihood(MascotPreparedInput prepared, double[] theta) {
        return evaluate(prepared, theta, null, false, false, false, false).logLikelihood;
    }

    public MascotLikelihoodDelegate.Result likelihoodAndGradient(MascotPreparedInput prepared, double[] theta) {
        return evaluate(prepared, theta, null, true, false, false, false);
    }

    public MascotLikelihoodDelegate.Result likelihoodAndAncestralStates(MascotPreparedInput prepared, double[] theta,
                                                                        double[] branchRates,
                                                                        boolean checkProbabilities) {
        return evaluate(prepared, theta, branchRates, false, true, checkProbabilities, false);
    }

    /**
     * Allocating wrapper for tests and callers that want a {@link Result}.
     * {@code branchRates} index biological lineage ids and scale only migration.
     * {@code copyFinalState} controls root-probability and active-lineage copies.
     */
    public MascotLikelihoodDelegate.Result evaluate(MascotPreparedInput prepared, double[] theta, double[] branchRates,
                                                    boolean computeGradient, boolean computeAncestralStates,
                                                    boolean checkProbabilities, boolean copyFinalState) {
        return evaluate(prepared, theta, branchRates, computeGradient, computeAncestralStates, checkProbabilities,
                copyFinalState, null);
    }

    /**
     * Test-only finite-difference hook: optional node/state log weights are applied
     * before coalescent-parent normalization. Package-private by design.
     */
    MascotLikelihoodDelegate.Result evaluateWithNodeLogWeightsForTesting(MascotPreparedInput prepared, double[] theta,
                                                                         double[] branchRates,
                                                                         double[] nodeLogWeights,
                                                                         boolean checkProbabilities) {
        return evaluate(prepared, theta, branchRates, false, false, checkProbabilities, false, nodeLogWeights);
    }

    private MascotLikelihoodDelegate.Result evaluate(MascotPreparedInput prepared, double[] theta,
                                                     double[] branchRates, boolean computeGradient,
                                                     boolean computeAncestralStates, boolean checkProbabilities,
                                                     boolean copyFinalState, double[] nodeLogWeights) {
        double[] gradientOut = computeGradient ? new double[parameterCount] : null;
        double[] clockGradientOut = computeGradient && branchRates != null
                ? new double[prepared.maxLineageId + 1] : null;
        double[] ancestralStateScoresOut = computeAncestralStates
                ? new double[(prepared.maxLineageId + 1) * stateCount] : null;

        double logLikelihood = evaluateCore(prepared, theta, branchRates, checkProbabilities, nodeLogWeights,
                gradientOut, clockGradientOut, ancestralStateScoresOut);

        MascotRuntime.ActiveState state = activeState;
        double[] rootProbabilities = copyFinalState ? state.copyProbabilities() : null;
        int[] activeLineages = copyFinalState ? state.copyActiveIds() : null;
        return new MascotLikelihoodDelegate.Result(logLikelihood, gradientOut, clockGradientOut,
                rootProbabilities, activeLineages, ancestralStateScoresOut);
    }

    /**
     * Caller-owned-output production path. Null outputs are skipped except that
     * gradient and ancestral-state requests share one reverse pass.
     */
    @Override
    public double calculateLikelihoodAndDerivatives(MascotPreparedInput prepared, double[] theta, double[] branchRates,
                                                    double[] gradientOut, double[] clockGradientOut,
                                                    double[] ancestralStateScoresOut, boolean checkProbabilities) {
        return evaluateCore(prepared, theta, branchRates, checkProbabilities, null,
                gradientOut, clockGradientOut, ancestralStateScoresOut);
    }

    /** Scalar-only evaluation: never builds the reverse tape, never allocates a {@link Result}. */
    @Override
    public double calculateLikelihood(MascotPreparedInput prepared, double[] theta, double[] branchRates,
                                      boolean checkProbabilities) {
        return evaluateCore(prepared, theta, branchRates, checkProbabilities, null, null, null, null);
    }

    private double evaluateCore(MascotPreparedInput prepared, double[] theta, double[] branchRates,
                                boolean checkProbabilities, double[] nodeLogWeights,
                                double[] gradientOut, double[] clockGradientOut, double[] ancestralStateScoresOut) {
        if (prepared == null) {
            throw new IllegalArgumentException("prepared input must not be null");
        }
        if (prepared.tipData.stateCount != stateCount) {
            throw new IllegalArgumentException("tip data stateCount " + prepared.tipData.stateCount +
                    " does not match GenericMascotLikelihoodDelegate stateCount " + stateCount);
        }
        if (theta == null || theta.length != parameterCount) {
            throw new IllegalArgumentException("theta dimension " + (theta == null ? -1 : theta.length) +
                    " does not match expected dimension " + parameterCount);
        }
        if (branchRates != null && branchRates.length <= prepared.maxLineageId) {
            throw new IllegalArgumentException("branchRates dimension " + branchRates.length +
                    " does not cover the maximum lineage id " + prepared.maxLineageId);
        }
        if (nodeLogWeights != null) {
            int expectedLength = (prepared.maxLineageId + 1) * stateCount;
            if (nodeLogWeights.length != expectedLength) {
                throw new IllegalArgumentException("nodeLogWeights dimension " + nodeLogWeights.length +
                        " does not match expected dimension " + expectedLength);
            }
        }
        if (gradientOut != null && gradientOut.length != parameterCount) {
            throw new IllegalArgumentException("gradientOut dimension " + gradientOut.length +
                    " does not match expected dimension " + parameterCount);
        }
        int expectedLineageDimension = prepared.maxLineageId + 1;
        if (clockGradientOut != null) {
            if (branchRates == null) {
                throw new IllegalArgumentException("clockGradientOut was requested but branchRates is null");
            }
            if (clockGradientOut.length != expectedLineageDimension) {
                throw new IllegalArgumentException("clockGradientOut dimension " + clockGradientOut.length +
                        " does not match expected dimension " + expectedLineageDimension);
            }
        }
        if (ancestralStateScoresOut != null && ancestralStateScoresOut.length != expectedLineageDimension * stateCount) {
            throw new IllegalArgumentException("ancestralStateScoresOut dimension " + ancestralStateScoresOut.length +
                    " does not match expected dimension " + (expectedLineageDimension * stateCount));
        }

        for (int epoch = 0; epoch < epochCount; epoch++) {
            updateEpochRates(theta, epoch, epochRates[epoch]);
        }

        MascotRuntime.ActiveState state = activeState;
        boolean runReverse = gradientOut != null || ancestralStateScoresOut != null;
        MascotAdjointTape.Store operations = null;
        if (runReverse) {
            operations = adjointEvaluator.resetTape(prepared.schedule.getIntervalCount() * 2 + epochCount + 1);
        }

        forwardEvaluator.forward(state, epochRates, workspace, prepared, branchRates, checkProbabilities,
                operations, nodeLogWeights);

        if (ancestralStateScoresOut != null) {
            Arrays.fill(ancestralStateScoresOut, 0, expectedLineageDimension * stateCount, Double.NaN);
        }

        if (runReverse) {
            adjointEvaluator.reverseInto(workspace, epochRates, operations, state.activeCount,
                    expectedLineageDimension, branchRates != null,
                    gradientOut, clockGradientOut, ancestralStateScoresOut);
        }

        return state.logLikelihood;
    }

    // ------------------------------------------------------------------
    // Epoch rate cache
    // ------------------------------------------------------------------

    private void updateEpochRates(double[] theta, int epoch, MascotRuntime.EpochRates rates) {
        int thetaOffset = epoch * parametersPerEpoch;
        int index = 0;
        for (int source = 0; source < stateCount; source++) {
            double rowSum = 0.0;
            int row = source * stateCount;
            for (int sink = 0; sink < stateCount; sink++) {
                if (source == sink) {
                    continue;
                }
                double rate = theta[thetaOffset + index];
                // Zero rates are valid (e.g. inactive BSSVS indicators);
                // only negative or non-finite rates break the equations.
                if (!(rate >= 0.0) || !Double.isFinite(rate)) {
                    throw new MascotLikelihoodDelegate.NumericalException("invalid migration rate for epoch " + epoch +
                            ", source " + source + ", sink " + sink + ": " + rate);
                }
                rates.migrationRates[index] = rate;
                rates.migrationMatrix[row + sink] = rate;
                rowSum += rate;
                if (!Double.isFinite(rowSum)) {
                    throw new MascotLikelihoodDelegate.NumericalException("invalid migration row sum for epoch " + epoch +
                            ", source " + source + ": " + rowSum);
                }
                index++;
            }
            rates.migrationMatrix[row + source] = -rowSum;
        }

        int etaOffset = thetaOffset + migrationParametersPerEpoch;
        for (int state = 0; state < stateCount; state++) {
            double logPopulation = theta[etaOffset + state];
            double inversePopulation = Math.exp(-logPopulation);
            if (!Double.isFinite(logPopulation) || !Double.isFinite(inversePopulation)) {
                throw new MascotLikelihoodDelegate.NumericalException("invalid log population size for epoch " + epoch +
                        ", state " + state + ": " + logPopulation);
            }
            rates.inversePopulation[state] = inversePopulation;
        }
    }

}
