/*
 * TimeTreeSearch.java
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

import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.Tree;
import dr.evolution.tree.treemetrics.SplitDistance;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.Reportable;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.logging.Logger;

/**
 * A deterministic topology search for TreeLime's time trees (plan 7), run outside MCMC: a move is accepted only if the
 * objective F (the maximizer's own, MAP or ML) rises by more than a threshold. F at any point reached on a candidate
 * topology is a lower bound on that topology's optimum, so accepting on it never lowers F; every optimizer call is
 * also guarded, and reverted if it ends lower.
 * <p>
 * This is the first version (M1a): a hill-climb by rooted NNIs with height repair. Each round scores every NNI
 * that keeps the constraints (T0: one evaluation after the move and the rate warm start; T1: the move's heights and
 * rates optimized locally when T0 is close), accepts the best one by replay, and re-fits the heights and rates (T3a,
 * the round maximizer). When no move gains, the top candidates get a short re-fit of all heights and rates (T2) as a
 * second look. A new best re-fits everything (T3b, the final maximizer). At the end the best state is restored and
 * re-fitted with the final maximizer's own settings.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class TimeTreeSearch implements Reportable {

    public enum StopReason {
        CONVERGED("converged: no NNI gains"),
        MAX_ROUNDS("maxRounds"),
        MAX_RUN_TIME("maxRunTime"),
        MAX_EVALUATIONS("maxEvaluations"),
        ERROR("error");

        private final String text;

        StopReason(String text) {
            this.text = text;
        }

        public String toString() {
            return text;
        }
    }

    private static final Logger LOGGER = Logger.getLogger("dr.evomodel");
    private static final double GUARD_TOLERANCE = 1E-6;

    private final TreeModel tree;
    private final Likelihood objective;
    private final MaximizerWrtParameter roundMaximizer;
    private final MaximizerWrtParameter finalMaximizer;
    private final ArbitraryBranchRates rates;
    private final TopologyConstraints constraints;
    private final SearchSettings settings;
    private final TrialEvaluator evaluator;
    private final List<Parameter> snapshotParameters;
    private final SearchLog log;

    // the latest run
    private boolean executed = false;
    private StopReason stopReason = null;
    private int rounds = 0;
    private int accepted = 0;
    private int guardReverts = 0;
    private long maximizerEvaluations = 0;
    private double fStart = Double.NaN;
    private double fRefitted = Double.NaN;
    private double fBest = Double.NaN;
    private double fFinal = Double.NaN;
    private long startNanos;
    private double seconds = 0.0;
    private double trialSeconds = 0.0;
    private double t3aSeconds = 0.0;
    private double t3bSeconds = 0.0;
    private Tree startTree = null;
    private Tree bestTree = null;

    public TimeTreeSearch(TreeModel tree, Likelihood objective, MaximizerWrtParameter roundMaximizer,
                          MaximizerWrtParameter finalMaximizer, ArbitraryBranchRates rates,
                          TopologyConstraints constraints, SearchSettings settings) {
        this.tree = tree;
        this.objective = objective;
        this.roundMaximizer = roundMaximizer;
        this.finalMaximizer = finalMaximizer;
        this.constraints = constraints;
        this.settings = settings;

        List<Parameter> roundParameters = new ArrayList<Parameter>();
        TrialEvaluator.flatten(roundMaximizer.getParameter(), roundParameters);
        // T1 and the warm start change branch rates: only when the round maximizer's parameters, which trials
        // record and replay, include them
        this.rates = rates != null && roundParameters.contains(rates.getRateParameter()) ? rates : null;
        BranchRateWarmStart warmStart = new BranchRateWarmStart(tree, this.rates, settings.rateWarmStart);
        LocalOptimizer local = settings.t1Sweeps > 0 ? new LocalOptimizer(tree, this.rates, settings.t1Sweeps) : null;
        this.evaluator = new TrialEvaluator(tree, objective, roundMaximizer, warmStart, local,
                settings.t2Iterations, settings.verifyEvery);

        this.snapshotParameters = new ArrayList<Parameter>();
        snapshotParameters.add(roundMaximizer.getParameter());
        if (finalMaximizer != null) {
            snapshotParameters.add(finalMaximizer.getParameter());
        }
        this.log = new SearchLog(tree, this.rates, objective);
    }

    public void run() {
        if (settings.seed >= 0) {
            MathUtils.setSeed(settings.seed);
        }
        executed = true;
        stopReason = null;
        rounds = 0;
        accepted = 0;
        guardReverts = 0;
        maximizerEvaluations = 0;
        trialSeconds = t3aSeconds = t3bSeconds = 0.0;
        log.clear();
        startNanos = System.nanoTime();
        startTree = new FlexibleTree(tree);

        if (!constraints.areSatisfied(tree)) {
            throw new IllegalStateException("The starting tree breaks a topology constraint");
        }

        TreeStateSnapshot best = null;
        try {
            fStart = evaluator.evaluate();
            log.row(0, "start", 0, 0, 0, 0, fStart, fStart, fStart, null, Double.NaN, Double.NaN, Double.NaN,
                    startTree, seconds(), evaluations(), "", guardReverts);
            say("start: F = " + format(fStart));

            double f = refit(finalMaximizer, startBudget(), "t3b", 0, fStart);
            fRefitted = f;
            fBest = f;
            best = snapshot();

            final int maxRounds = settings.maxRounds > 0 ? settings.maxRounds : tree.getInternalNodeCount();
            stopReason = StopReason.MAX_ROUNDS;
            for (int round = 1; round <= maxRounds; ++round) {
                StopReason budget = checkBudgets();
                if (budget != null) {
                    stopReason = budget;
                    break;
                }
                rounds = round;
                final double before = f;
                final boolean changed = nniRound(round);
                if (!changed) {
                    stopReason = checkBudgets() != null ? checkBudgets() : StopReason.CONVERGED;
                    break;
                }
                accepted++;
                f = refit(roundMaximizer, t3aBudget(), "t3a", round, evaluator.evaluate());
                if (f > fBest + settings.modelEpsilon) {
                    f = refit(finalMaximizer, bestBudget(), "t3b", round, f);
                }
                if (f > fBest) {
                    fBest = f;
                    best = snapshot();
                }
                if (f - before < settings.epsilon) {
                    stopReason = StopReason.CONVERGED;
                    break;
                }
            }

            finish(best);
        } catch (RuntimeException e) {
            stopReason = StopReason.ERROR;
            if (best != null) {
                try {
                    best.restore();
                } catch (RuntimeException ignored) {
                    // report the original error
                }
            }
            throw e;
        } finally {
            seconds = seconds();
        }
    }

    private void finish(TreeStateSnapshot best) {
        if (evaluator.evaluate() < fBest) {
            best.restore();
        }
        // the final fit, with the final maximizer's own settings
        fFinal = refit(finalMaximizer, null, "t3b_final", rounds, evaluator.evaluate());
        fBest = Math.max(fBest, fFinal);
        bestTree = new FlexibleTree(tree);
        seconds = seconds();
        say("stopped (" + stopReason + ") after " + rounds + " rounds: F = " + format(fFinal) + " (" +
                format(fFinal - fStart) + " over the start), " + accepted + " moves accepted, " +
                String.format("%.1fs", seconds));
        try {
            if (settings.logFileName != null) {
                log.writeLog(new File(settings.logFileName));
            }
            if (settings.treeFileName != null) {
                log.writeTree(new File(settings.treeFileName), fFinal, stopReason.toString());
            }
        } catch (IOException e) {
            throw new RuntimeException("Could not write the search output: " + e.getMessage(), e);
        }
    }

    /**
     * One round: score every allowed NNI, accept the best if it gains, else give the closest candidates a T2 look.
     *
     * @return whether a move was accepted
     */
    private boolean nniRound(int round) {
        final long started = System.nanoTime();
        final double fInc = evaluator.evaluate();
        final long t1Before = evaluator.getT1Count();
        final long t2Before = evaluator.getT2Count();
        final List<RootedNniMove> moves = NniNeighbourhood.enumerate(tree, constraints);

        List<TrialEvaluator.Trial> trials = new ArrayList<TrialEvaluator.Trial>();
        for (RootedNniMove move : moves) {
            if (checkBudgets() != null) {
                break;
            }
            trials.add(evaluator.trial(move, fInc - settings.deltaT1, false));
        }
        TrialEvaluator.Trial winner = null;
        for (TrialEvaluator.Trial trial : trials) {
            if (trial.best > fInc + settings.moveEpsilon && (winner == null || trial.best > winner.best)) {
                winner = trial;
            }
        }

        if (winner == null && checkBudgets() == null) {
            // a second look at the closest candidates with all heights and rates re-fitted (T2)
            List<TrialEvaluator.Trial> close = new ArrayList<TrialEvaluator.Trial>();
            for (TrialEvaluator.Trial trial : trials) {
                if (trial.best > fInc - settings.deltaKeep) {
                    close.add(trial);
                }
            }
            Collections.sort(close, new Comparator<TrialEvaluator.Trial>() {
                public int compare(TrialEvaluator.Trial x, TrialEvaluator.Trial y) {
                    return Double.compare(y.best, x.best); // stable: ties keep the enumeration order
                }
            });
            for (int k = 0; k < Math.min(settings.t2Candidates, close.size()); ++k) {
                TrialEvaluator.Trial trial = evaluator.trial(close.get(k).move, Double.NEGATIVE_INFINITY, true);
                if (trial.best > fInc + settings.moveEpsilon && (winner == null || trial.best > winner.best)) {
                    winner = trial;
                }
            }
        }

        final int nT0 = trials.size();
        final int nT1 = (int) (evaluator.getT1Count() - t1Before);
        final int nT2 = (int) (evaluator.getT2Count() - t2Before);
        double f = fInc;
        if (winner != null) {
            f = evaluator.accept(winner);
        }
        trialSeconds += (System.nanoTime() - started) / 1E9;
        log.row(round, "nni", nT0, nT1, nT2, winner != null ? 1 : 0, fInc, f, Math.max(fBest, f),
                winner != null ? winner.move.describe() : null,
                winner != null ? winner.f0 - fInc : Double.NaN,
                winner != null ? winner.f1 - fInc : Double.NaN,
                winner != null ? winner.f2 - fInc : Double.NaN,
                startTree, seconds(), evaluations(), "", guardReverts);
        say("round " + round + " NNI: F = " + format(f) + " (" + format(f - fInc) + ") " +
                (winner != null ? winner.move.describe() : "no move gains") + ", " + nT0 + " moves, " + nT1 +
                " T1, " + nT2 + " T2, " + String.format("%.1fs", seconds()));
        return winner != null;
    }

    // a maximizer call, reverted if it ends below where it started
    private double refit(MaximizerWrtParameter maximizer, MaximizerWrtParameter.Budget budget, String phase,
                         int round, double before) {
        if (maximizer == null) {
            return before;
        }
        final long started = System.nanoTime();
        final TreeStateSnapshot snapshot = snapshot();
        maximizer.maximize(budget);
        maximizerEvaluations += maximizer.getValueCalls();
        double f = evaluator.evaluate();
        if (f < before - GUARD_TOLERANCE * Math.max(1.0, Math.abs(before))) {
            snapshot.restore();
            guardReverts++;
            f = evaluator.evaluate();
        }
        final double elapsed = (System.nanoTime() - started) / 1E9;
        if (phase.startsWith("t3a")) {
            t3aSeconds += elapsed;
        } else {
            t3bSeconds += elapsed;
        }
        log.row(round, phase, 0, 0, 0, 0, before, f, Double.isNaN(fBest) ? f : Math.max(fBest, f), null, Double.NaN,
                Double.NaN, Double.NaN, startTree, seconds(), evaluations(), String.valueOf(maximizer.getReason()), guardReverts);
        return f;
    }

    private TreeStateSnapshot snapshot() {
        return new TreeStateSnapshot(tree, snapshotParameters, objective);
    }

    private MaximizerWrtParameter.Budget startBudget() {
        return new MaximizerWrtParameter.Budget(settings.t3bStartIterations)
                .setLimits(0.0, 0, settings.t3bStartRelTol, 10).setQuiet(true);
    }

    private MaximizerWrtParameter.Budget bestBudget() {
        return new MaximizerWrtParameter.Budget(settings.t3bBestIterations)
                .setLimits(0.0, 0, settings.t3bBestRelTol, 10).setQuiet(true);
    }

    private MaximizerWrtParameter.Budget t3aBudget() {
        return new MaximizerWrtParameter.Budget(settings.t3aIterations)
                .setLimits(0.0, 0, settings.t3aRelTol, settings.t3aStallWindow).setQuiet(true);
    }

    private StopReason checkBudgets() {
        if (settings.maxRunTime > 0 && seconds() >= settings.maxRunTime) {
            return StopReason.MAX_RUN_TIME;
        }
        if (settings.maxEvaluations > 0 && evaluations() >= settings.maxEvaluations) {
            return StopReason.MAX_EVALUATIONS;
        }
        return null;
    }

    private long evaluations() {
        return evaluator.getEvaluations() + maximizerEvaluations;
    }

    private double seconds() {
        return (System.nanoTime() - startNanos) / 1E9;
    }

    private void say(String message) {
        if (settings.printToScreen) {
            LOGGER.info("[search] " + message);
        }
    }

    private static String format(double value) {
        return String.format("%.4f", value);
    }

    // ---------------------------------------------------------------------------------------------------------------

    public boolean wasExecuted() {
        return executed;
    }

    public StopReason getStopReason() {
        return stopReason;
    }

    public int getRounds() {
        return rounds;
    }

    public int getAcceptedMoves() {
        return accepted;
    }

    public int getGuardReverts() {
        return guardReverts;
    }

    public double getStartValue() {
        return fStart;
    }

    public double getBestValue() {
        return fBest;
    }

    public double getFinalValue() {
        return fFinal;
    }

    public TrialEvaluator getEvaluator() {
        return evaluator;
    }

    public SearchLog getLog() {
        return log;
    }

    public TopologyConstraints getConstraints() {
        return constraints;
    }

    public TreeModel getTree() {
        return tree;
    }

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder("Topology search (rooted NNI hill-climb)\n");
        if (!executed) {
            return sb.append("Not yet executed.\n").toString();
        }
        sb.append("Settings: ").append(settings.describe()).append("\n");
        sb.append("StopReason: ").append(stopReason).append("\n");
        sb.append("Rounds: ").append(rounds).append("\n");
        sb.append("AcceptedMoves: ").append(accepted).append("\n");
        sb.append("Trials: ").append(evaluator.getTrials()).append(" (T1 ").append(evaluator.getT1Count())
                .append(", T2 ").append(evaluator.getT2Count()).append(")\n");
        sb.append("GuardReverts: ").append(guardReverts).append("\n");
        if (settings.verifyEvery > 0) {
            sb.append("Verifications: ").append(evaluator.getVerifications()).append("\n");
        }
        sb.append("Evaluations: ").append(evaluations()).append("\n");
        sb.append("F: start ").append(fStart).append(", after the starting re-fit ").append(fRefitted)
                .append(", final ").append(fFinal).append(" (").append(fFinal - fStart).append(" over the start)\n");
        if (startTree != null && bestTree != null) {
            sb.append("RF start -> final: rooted ").append(SplitDistance.rooted(startTree, bestTree))
                    .append(", unrooted ").append(SplitDistance.unrooted(startTree, bestTree)).append("\n");
            sb.append("Topology: ").append(TopologyHash.toHex(TopologyHash.of(startTree))).append(" -> ")
                    .append(TopologyHash.toHex(TopologyHash.of(bestTree))).append("\n");
        }
        sb.append("Constraints: ").append(constraints.size());
        if (constraints.getUnfilteredCount() > 0) {
            sb.append(" (").append(constraints.getUnfilteredCount()).append(" more scored through F only)");
        }
        sb.append("\n");
        sb.append(String.format("Seconds: %.3f (trials %.3f, T3a %.3f, T3b %.3f)%n", seconds, trialSeconds,
                t3aSeconds, t3bSeconds));
        if (settings.logFileName != null) {
            sb.append("Log: ").append(settings.logFileName).append("\n");
        }
        if (settings.treeFileName != null) {
            sb.append("Tree: ").append(settings.treeFileName).append("\n");
        }
        return sb.toString();
    }
}
