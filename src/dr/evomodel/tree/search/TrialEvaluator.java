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
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.evomodel.treedatalikelihood.discrete.NodeHeightProxyParameter;
import dr.inference.model.AbstractModel;
import dr.inference.model.CompoundModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;

import java.util.ArrayList;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

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
    public static final class Theta {
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
    private final List<BeagleDataLikelihoodDelegate> beagles = new ArrayList<BeagleDataLikelihoodDelegate>();

    private long evaluations = 0;
    private long trials = 0;
    private long t1Count = 0;
    private long t2Count = 0;
    private long verifications = 0;
    private long rescues = 0;

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
        for (Model m : models(model)) {
            if (m instanceof TreeDataLikelihood &&
                    ((TreeDataLikelihood) m).getDataLikelihoodDelegate() instanceof BeagleDataLikelihoodDelegate) {
                beagles.add((BeagleDataLikelihoodDelegate) ((TreeDataLikelihood) m).getDataLikelihoodDelegate());
            }
        }
    }

    /**
     * @return every model reachable from the root
     */
    public static Set<Model> models(Model root) {
        Map<Model, Boolean> seen = new IdentityHashMap<Model, Boolean>();
        List<Model> stack = new ArrayList<Model>();
        stack.add(root);
        while (!stack.isEmpty()) {
            Model model = stack.remove(stack.size() - 1);
            if (model == null || seen.containsKey(model)) {
                continue;
            }
            seen.put(model, Boolean.TRUE);
            if (model instanceof AbstractModel) {
                AbstractModel abstractModel = (AbstractModel) model;
                for (int i = 0; i < abstractModel.getModelCount(); ++i) {
                    stack.add(abstractModel.getModel(i));
                }
            } else if (model instanceof CompoundModel) {
                CompoundModel compound = (CompoundModel) model;
                for (int i = 0; i < compound.getModelCount(); ++i) {
                    stack.add(compound.getModel(i));
                }
            }
        }
        return seen.keySet();
    }

    // back to the stored state. If BEAGLE rescaled inside the window it overwrote stored partials (and, the first
    // time, left no stored scale factors), so everything is recomputed now, before the next window stores this state
    private void restore(long underflowsBefore) {
        model.restoreModelState();
        if (underflows() != underflowsBefore) {
            ++rescues;
            objective.makeDirty();
            evaluate();
        }
    }

    private long underflows() {
        long count = 0;
        for (BeagleDataLikelihoodDelegate beagle : beagles) {
            count += beagle.getUnderflowCount();
        }
        return count;
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
        return trial(move, t1Threshold, t2 ? t2Budget : null);
    }

    /**
     * Scores a move as above, with T2 for the given budget (null: no T2).
     */
    public Trial trial(TimeTreeMove move, double t1Threshold, MaximizerWrtParameter.Budget t2Budget) {
        final double before = evaluate(); // TreeDataLikelihood.storeState needs a computed likelihood
        final long underflows = underflows();
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
            if (t2Budget != null && roundMaximizer != null) {
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
            restore(underflows);
        }
        ++trials;
        if (verifyEvery > 0 && trials % verifyEvery == 0) {
            verify(before);
        }
        return trial;
    }

    /**
     * An edit kept only if it says so and leaves a finite objective (perturbation by BEAST's operators).
     */
    public interface Edit {
        /**
         * Changes the tree (in its own edit transactions).
         *
         * @return false to reject the edit outright
         */
        boolean apply();
    }

    /**
     * Applies an edit inside a store window and keeps it if it succeeds and the objective stays finite, else restores
     * the state as it was (recomputing everything if BEAGLE rescaled meanwhile).
     *
     * @return whether the edit was kept
     */
    public boolean tryEdit(Edit edit) {
        evaluate();
        final long underflows = underflows();
        model.storeModelState();
        boolean kept = false;
        try {
            if (edit.apply()) {
                final double f = evaluate();
                kept = f > Double.NEGATIVE_INFINITY && !Double.isNaN(f);
            }
        } finally {
            if (kept) {
                model.acceptModelState();
            } else {
                restore(underflows);
            }
        }
        return kept;
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
        final long underflows = underflows();
        model.storeModelState();
        tree.beginTreeEdit();
        trial.move.applyTopology(tree);
        setHeights(trial.theta.heights);
        tree.endTreeEdit();
        setValues(trial.theta.values);
        final double f = evaluate();
        if (!(Math.abs(f - trial.best) <= REPLAY_TOLERANCE * Math.max(1.0, Math.abs(trial.best)))) {
            restore(underflows);
            throw new IllegalStateException("Replaying " + trial.move.describe() + " gives F = " + f +
                    ", not the trial's " + trial.best);
        }
        model.acceptModelState();
        return f;
    }

    /**
     * Accepts several scored moves at once (plan 7 §7.3 batch). The moves must not touch each other's quartets: each
     * one changed only its own local heights and rates from the base state, so their changes, applied together on
     * the base state, do not overlap.
     *
     * @param base the state every trial started from
     * @return F after the batch, or NaN (with the state as it was) when the merged heights are not a valid time tree
     */
    public double acceptBatch(List<Trial> batch, Theta base) {
        final double[] heights = base.heights.clone();
        final double[][] values = new double[base.values.length][];
        for (int p = 0; p < values.length; ++p) {
            values[p] = base.values[p].clone();
        }
        for (Trial trial : batch) {
            for (int i = 0; i < heights.length; ++i) {
                if (!Double.isNaN(heights[i]) && trial.theta.heights[i] != base.heights[i]) {
                    heights[i] = trial.theta.heights[i];
                }
            }
            for (int p = 0; p < values.length; ++p) {
                for (int i = 0; i < values[p].length; ++i) {
                    if (trial.theta.values[p][i] != base.values[p][i]) {
                        values[p][i] = trial.theta.values[p][i];
                    }
                }
            }
        }

        evaluate();
        final long underflows = underflows();
        model.storeModelState();
        tree.beginTreeEdit();
        for (Trial trial : batch) {
            if (!trial.move.isValid(tree)) {
                tree.endTreeEdit();
                restore(underflows);
                throw new IllegalStateException("The batch's moves overlap at " + trial.move.describe());
            }
            trial.move.applyTopology(tree);
        }
        setHeights(heights);
        tree.endTreeEdit();
        setValues(values);
        final double f = isTimeTree() ? evaluate() : Double.NaN;
        if (Double.isNaN(f) || Double.isInfinite(f)) {
            restore(underflows);
            return Double.NaN;
        }
        model.acceptModelState();
        return f;
    }

    private boolean isTimeTree() {
        for (int i = 0; i < tree.getNodeCount(); ++i) {
            final NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node) && !(tree.getNodeHeight(node) < tree.getNodeHeight(tree.getParent(node)))) {
                return false;
            }
        }
        return true;
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

    /**
     * @return the current heights and round-parameter values
     */
    public Theta capture() {
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

    public MaximizerWrtParameter.Budget getT2Budget() {
        return t2Budget;
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

    /**
     * @return the trials after which the state was recomputed because BEAGLE rescaled during the trial
     */
    public long getUnderflowRecomputes() {
        return rescues;
    }

    public long getVerifications() {
        return verifications;
    }

    public List<Parameter> getRoundParameters() {
        return roundParameters;
    }
}
