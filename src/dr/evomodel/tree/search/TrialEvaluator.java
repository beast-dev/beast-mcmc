/*
 * TrialEvaluator.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodel.tree.search;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.evomodel.treedatalikelihood.discrete.NodeHeightProxyParameter;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;

import java.util.ArrayList;
import java.util.List;

/**
 * Scores topology moves inside TreeModel's own store window, and accepts one by replaying it (plan 7 §4.4.1):
 * <pre>
 *     store; apply the move (repair, rate warm start); T0 = F; [T1: local heights and rates]; [T2: round maximizer];
 *     record the move with its best F and the heights and round parameters there; restore
 * </pre>
 * Only the tree's heights and the round maximizer's parameters (rates, a strict clock's rate) change inside the
 * window, so the restore brings everything back; BEAGLE only swaps buffer indices. To accept, the move is applied again
 * with the recorded values and F must come back to the recorded value. With verifyEvery N, every N-th trial is
 * followed by a full recomputation that must give the pre-trial F.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public final class TrialEvaluator {

    /**
     * A scored move and the state (heights by node number, round-parameter values) at its best F.
     */
    public static final class Trial {
        public final TimeTreeMove move;
        public final double f0;
        public double f1 = Double.NaN;
        public double f2 = Double.NaN;
        public double best;
        Theta theta;

        Trial(TimeTreeMove move, double f0) {
            this.move = move;
            this.f0 = f0;
            this.best = f0;
        }

        /**
         * @return the tier of the best value: 0, 1 or 2
         */
        public int getTier() {
            return best == f2 ? 2 : (best == f1 ? 1 : 0);
        }
    }

    // the values the trials change: internal heights by node number, and the round parameters' values
    static final class Theta {
        final double[] heights;
        final double[][] values;

        Theta(double[] heights, double[][] values) {
            this.heights = heights;
            this.values = values;
        }
    }

    private static final double REPLAY_TOLERANCE = 1E-6;
    private static final double VERIFY_TOLERANCE = 1E-8;

    private final TreeModel tree;
    private final Likelihood objective;
    private final Model model;
    private final MaximizerWrtParameter roundMaximizer;
    private final List<Parameter> roundParameters; // without node-height proxies
    private final BranchRateWarmStart warmStart;
    private final LocalOptimizer localOptimizer;
    private final int verifyEvery;
    private final MaximizerWrtParameter.Budget t2Budget;

    private long evaluations = 0;
    private long trials = 0;
    private long t1Count = 0;
    private long t2Count = 0;
    private long verifications = 0;

    public TrialEvaluator(TreeModel tree, Likelihood objective, MaximizerWrtParameter roundMaximizer,
                          BranchRateWarmStart warmStart, LocalOptimizer localOptimizer, int t2Iterations,
                          int verifyEvery) {
        this.tree = tree;
        this.objective = objective;
        this.model = objective.getModel();
        this.roundMaximizer = roundMaximizer;
        this.roundParameters = new ArrayList<Parameter>();
        if (roundMaximizer != null) {
            flatten(roundMaximizer.getParameter(), roundParameters);
        }
        this.warmStart = warmStart;
        this.localOptimizer = localOptimizer;
        this.verifyEvery = verifyEvery;
        this.t2Budget = new MaximizerWrtParameter.Budget(t2Iterations).setQuiet(true);
    }

    static void flatten(Parameter parameter, List<Parameter> leaves) {
        if (parameter instanceof CompoundParameter) {
            CompoundParameter compound = (CompoundParameter) parameter;
            for (int i = 0; i < compound.getParameterCount(); ++i) {
                flatten(compound.getParameter(i), leaves);
            }
        } else if (!(parameter instanceof NodeHeightProxyParameter)) {
            leaves.add(parameter);
        }
    }

    /**
     * @return the objective at the current state
     */
    public double evaluate() {
        ++evaluations;
        return objective.getLogLikelihood();
    }

    /**
     * Scores a move: T0 always, T1 when T0 exceeds t1Threshold (and a local optimizer is set), T2 when asked. The
     * tree and every parameter are back as they were afterwards.
     */
    public Trial trial(TimeTreeMove move, double t1Threshold, boolean t2) {
        final double before = evaluate(); // TreeDataLikelihood.storeState needs a computed likelihood
        model.storeModelState();
        final Trial trial;
        try {
            move.apply(tree, warmStart);
            trial = new Trial(move, evaluate());
            trial.theta = capture();
            if (localOptimizer != null && trial.f0 > t1Threshold) {
                ++t1Count;
                trial.f1 = localOptimizer.optimize(new LocalOptimizer.Objective() {
                    public double value() {
                        return evaluate();
                    }
                }, trial.f0, move.localHeights(tree), move.localRates(tree));
                if (trial.f1 > trial.best) {
                    trial.best = trial.f1;
                    trial.theta = capture();
                }
            }
            if (t2 && roundMaximizer != null) {
                ++t2Count;
                // from the T1 point as it is in the tree, with a cold L-BFGS memory
                roundMaximizer.maximize(t2Budget);
                trial.f2 = evaluate();
                if (trial.f2 > trial.best) {
                    trial.best = trial.f2;
                    trial.theta = capture();
                }
            }
        } finally {
            model.restoreModelState();
        }
        ++trials;
        if (verifyEvery > 0 && trials % verifyEvery == 0) {
            verify(before);
        }
        return trial;
    }

    /**
     * Accepts a scored move: applies it again with the recorded heights and parameters, checks that F is the
     * recorded one, and keeps the state.
     *
     * @return F after the replay
     * @throws IllegalStateException when the replay gives another F (the state is restored first)
     */
    public double accept(Trial trial) {
        evaluate();
        model.storeModelState();
        tree.beginTreeEdit();
        trial.move.applyTopology(tree);
        setHeights(trial.theta.heights);
        tree.endTreeEdit();
        setValues(trial.theta.values);
        final double f = evaluate();
        if (!(Math.abs(f - trial.best) <= REPLAY_TOLERANCE * Math.max(1.0, Math.abs(trial.best)))) {
            model.restoreModelState();
            throw new IllegalStateException("Replaying " + trial.move.describe() + " gives F = " + f +
                    ", not the trial's " + trial.best);
        }
        model.acceptModelState();
        return f;
    }

    /**
     * Makes every model recompute from scratch and checks the result against the expected F.
     */
    public void verify(double expected) {
        ++verifications;
        objective.makeDirty();
        final double f = evaluate();
        if (!(Math.abs(f - expected) <= VERIFY_TOLERANCE * Math.max(1.0, Math.abs(expected)))) {
            throw new IllegalStateException("A full recomputation gives F = " + f + ", not " + expected +
                    ": an incremental update was missed");
        }
    }

    Theta capture() {
        final double[] heights = new double[tree.getNodeCount()];
        for (int i = 0; i < heights.length; ++i) {
            final NodeRef node = tree.getNode(i);
            heights[i] = tree.isExternal(node) ? Double.NaN : tree.getNodeHeight(node);
        }
        final double[][] values = new double[roundParameters.size()][];
        for (int p = 0; p < values.length; ++p) {
            values[p] = roundParameters.get(p).getParameterValues();
        }
        return new Theta(heights, values);
    }

    // inside the caller's edit transaction: the height events are queued with the topology's
    private void setHeights(double[] heights) {
        for (int i = 0; i < heights.length; ++i) {
            final NodeRef node = tree.getNode(i);
            if (!Double.isNaN(heights[i]) && tree.getNodeHeight(node) != heights[i]) {
                tree.setNodeHeight(node, heights[i]);
            }
        }
    }

    private void setValues(double[][] values) {
        for (int p = 0; p < values.length; ++p) {
            final Parameter parameter = roundParameters.get(p);
            for (int i = 0; i < values[p].length; ++i) {
                if (parameter.getParameterValue(i) != values[p][i]) {
                    parameter.setParameterValue(i, values[p][i]);
                }
            }
        }
    }

    public long getEvaluations() {
        return evaluations;
    }

    public long getTrials() {
        return trials;
    }

    public long getT1Count() {
        return t1Count;
    }

    public long getT2Count() {
        return t2Count;
    }

    public long getVerifications() {
        return verifications;
    }

    public List<Parameter> getRoundParameters() {
        return roundParameters;
    }
}
