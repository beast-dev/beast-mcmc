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
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.treemetrics.SplitDistance;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.Reportable;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

/**
 * A deterministic topology search for TreeLime's time trees (plan 7), run outside MCMC: a move is accepted only if the
 * objective F (the maximizer's own, MAP or ML) rises by more than a threshold. F at any point reached on a candidate
 * topology is a lower bound on that topology's optimum, so accepting on it never lowers F; every optimizer call is
 * also guarded, and reverted if it ends lower.
 * <p>
 * The search is a hill-climb by rooted NNIs with height repair. Each round scores the NNIs that keep the
 * constraints (T0: one evaluation after the move; T1: the move's heights and rates optimized locally when T0 is
 * close), only those near the last change when speedNNI is on (and all of them before stopping). It applies every
 * gaining move that does not touch another's neighbourhood at once, falling back to the best single move if that
 * batch ends lower, then re-fits the heights and rates (T3a, the round maximizer). When no move gains, the closest
 * candidates get a short re-fit of all heights and rates (T2), and near-ties a longer one. A new best re-fits
 * everything (T3b, the final maximizer) and tries the root positions one and two edges away. At the end the best
 * state is restored and re-fitted with the final maximizer's own settings.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class TimeTreeSearch implements Reportable, Citable {

    public enum StopReason {
        CONVERGED("converged: no NNI gains"),
        MAX_ROUNDS("maxRounds"),
        N_STOP("nStop: no new best in nStop perturbed hill-climbs"),
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
    private int secondChances = 0;
    private int rootScans = 0;
    private int rootMoves = 0;
    private int[] lastApplied = null; // the nodes of the last round's moves, for speedNNI
    private TreeStateSnapshot bestState = null;
    private final CandidateSet candidates;
    private final OperatorPerturbation perturbation; // null when nStop is 0
    private int climb = 0; // the hill-climb in progress: 0 from the starting tree, k after the k-th perturbation
    private int perturbations = 0;
    private double[] nniDelta = null;
    private int nniDeltaNegative = 0;
    private int sprRounds = 0;
    private int swapRounds = 0;
    private int swapAccepted = 0;
    private int rootSprRounds = 0;
    private int rootSprAccepted = 0;
    private int sprAccepted = 0;
    private double sprCutoff = -1.0; // the SPR scan's cutoff below the incumbent; < 0 before the first round
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
        this.candidates = new CandidateSet(settings.candidateSetSize);
        this.perturbation = settings.nStop > 0 ? new OperatorPerturbation(tree, evaluator, constraints,
                settings.perturbOperators, settings.leapSize) : null;
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
        secondChances = rootScans = rootMoves = sprRounds = sprAccepted = swapRounds = swapAccepted = 0;
        rootSprRounds = rootSprAccepted = 0;
        sprCutoff = -1.0;
        maximizerEvaluations = 0;
        trialSeconds = t3aSeconds = t3bSeconds = 0.0;
        log.clear();
        startNanos = System.nanoTime();
        startTree = new FlexibleTree(tree);

        if (!constraints.areSatisfied(tree)) {
            throw new IllegalStateException("The starting tree breaks a topology constraint");
        }
        if (settings.njStart && !settings.currentStart) {
            njStart(true);
        }

        bestState = null;
        candidates.clear();
        climb = 0;
        perturbations = 0;
        log.setClimb(0);
        try {
            fStart = evaluator.evaluate();
            log.row(0, "start", 0, 0, 0, 0, fStart, fStart, fStart, null, Double.NaN, Double.NaN, Double.NaN,
                    startTree, seconds(), evaluations(), "", guardReverts);
            say("start: F = " + format(fStart));
            // the second start begins from the parameters as they were, like a first one: after the first climb's
            // fit (a small clock scale, say) the NJ tree's rates can make F -Infinity and the re-fit stall far below
            final TreeStateSnapshot initialState = settings.njStart && settings.currentStart ? snapshot() : null;

            double f = refit(finalMaximizer, startBudget(), "t3b", 0, fStart);
            fRefitted = f;
            fBest = f;
            bestState = snapshot();

            stopReason = hillClimb(f, true);
            addCandidate();

            // the second start (§7.1): the NJ tree, fitted and climbed like the first; the perturbations then draw
            // from the best trees of both
            if (settings.njStart && settings.currentStart && checkBudgets() == null &&
                    (stopReason == StopReason.CONVERGED || stopReason == StopReason.MAX_ROUNDS)) {
                climb++;
                log.setClimb(climb);
                initialState.restore();
                if (njStart(false)) {
                    final double fNj = evaluator.evaluate();
                    log.row(rounds, "start", 0, 0, 0, 0, fNj, fNj, fBest, null, Double.NaN, Double.NaN, Double.NaN,
                            startTree, seconds(), evaluations(), "nj", guardReverts);
                    stopReason = hillClimb(refit(finalMaximizer, startBudget(), "t3b", rounds, fNj), true);
                    addCandidate();
                } else {
                    bestState.restore();
                }
            }

            // the perturbation phase (§7.2): a perturbed copy of a good tree, re-fitted and hill-climbed, until nStop
            // climbs in a row bring no new best
            if (perturbation != null &&
                    (stopReason == StopReason.CONVERGED || stopReason == StopReason.MAX_ROUNDS)) {
                final int moves = Math.max(1, (int) Math.floor(settings.perturbStrength *
                        (tree.getExternalNodeCount() - 2)));
                int unsuccessful = 0;
                while (unsuccessful < settings.nStop) {
                    final StopReason budget = checkBudgets();
                    if (budget != null) {
                        stopReason = budget;
                        break;
                    }
                    climb++;
                    log.setClimb(climb);
                    final double before = fBest;
                    final CandidateSet.Candidate parent = candidates.drawParent(settings.parentPool);
                    parent.getSnapshot().restore();
                    final double fParent = evaluator.evaluate();
                    final int kept = perturbation.perturb(moves);
                    perturbations++;
                    final double fPerturbed = evaluator.evaluate();
                    log.row(rounds, "perturb", kept, 0, 0, kept, fParent, fPerturbed, fBest, null, Double.NaN,
                            Double.NaN, Double.NaN, startTree, seconds(), evaluations(), "parent climb " +
                                    parent.climb + ", " + moves + " moves", guardReverts);
                    say("perturbation " + climb + ": F = " + format(fParent) + " -> " + format(fPerturbed) + ", " +
                            kept + " operator moves");
                    final double fitted = refit(roundMaximizer, t3aBudget(), "t3a", rounds, fPerturbed);
                    final StopReason climbed = hillClimb(fitted, false);
                    addCandidate();
                    if (climbed != StopReason.CONVERGED && climbed != StopReason.MAX_ROUNDS) {
                        stopReason = climbed;
                        break;
                    }
                    unsuccessful = fBest > before + settings.epsilon ? 0 : unsuccessful + 1;
                    if (unsuccessful >= settings.nStop) {
                        stopReason = StopReason.N_STOP;
                    }
                }
            }

            finish(bestState);
        } catch (RuntimeException e) {
            stopReason = StopReason.ERROR;
            if (bestState != null) {
                try {
                    bestState.restore();
                } catch (RuntimeException ignored) {
                    // report the original error
                }
            }
            throw e;
        } finally {
            seconds = seconds();
        }
    }

    /**
     * The hill-climb from the current state (§7.3): NNI rounds until they stop gaining, then SPR and subtree swaps
     * when on, then NNIs again; each change re-fits heights and rates (T3a), and a new best re-fits everything (T3b)
     * and scans the root positions nearby.
     *
     * @param f     the objective at the current state
     * @param start a start's climb, whose own best triggers the T3b re-fits, so that a second start below the
     *              incumbent climbs as it would alone; otherwise (a perturbed tree) only a new global best does
     * @return CONVERGED, MAX_ROUNDS (this climb's), or the budget that ran out
     */
    private StopReason hillClimb(double f, boolean start) {
        final int maxRounds = settings.maxRounds > 0 ? settings.maxRounds : tree.getInternalNodeCount();
        double fClimb = f; // this climb's best
        lastApplied = null;
        boolean beyondNniNext = false; // the NNI rounds have converged: SPR and swaps next
        for (int r = 1; r <= maxRounds; ++r) {
            StopReason budget = checkBudgets();
            if (budget != null) {
                return budget;
            }
            final int round = ++rounds;
            final double before = f;
            // NNIs until they stop gaining, then SPR, then subtree swaps, then NNIs again from the full tree
            boolean beyondNni = beyondNniNext;
            int applied = beyondNniNext ? 0 : nniRound(round);
            if (applied == 0 && settings.spr && checkBudgets() == null) {
                applied = sprPhase(round);
                lastApplied = null;
                beyondNni = true;
            }
            if (applied == 0 && settings.swap && checkBudgets() == null) {
                applied = swapRound(round);
                lastApplied = null;
                beyondNni = true;
            }
            if (applied == 0 && settings.rootMoves && checkBudgets() == null) {
                applied = rootSprRound(round);
                lastApplied = null;
                beyondNni = true;
            }
            beyondNniNext = false;
            if (applied == 0) {
                return checkBudgets() != null ? checkBudgets() : StopReason.CONVERGED;
            }
            accepted += applied;
            f = refit(roundMaximizer, t3aBudget(), "t3a", round, evaluator.evaluate());
            // a new best re-fits everything, then looks at the root positions nearby (§6.3)
            while (f > (start ? fClimb : fBest) + settings.modelEpsilon) {
                f = refit(finalMaximizer, bestBudget(), "t3b", round, f);
                fClimb = Math.max(fClimb, f);
                if (f > fBest) {
                    fBest = f;
                    bestState = snapshot();
                }
                if (!settings.rootScan || checkBudgets() != null || !rootScan(round)) {
                    break;
                }
                accepted++;
                f = refit(roundMaximizer, t3aBudget(), "t3a", round, evaluator.evaluate());
            }
            fClimb = Math.max(fClimb, f);
            if (f > fBest) {
                fBest = f;
                bestState = snapshot();
            }
            if (f - before < settings.epsilon) {
                if ((settings.spr || settings.swap || settings.rootMoves) && !beyondNni && checkBudgets() == null) {
                    beyondNniNext = true; // the NNIs gain too little to go on: the wider moves before stopping
                    continue;
                }
                return StopReason.CONVERGED;
            }
        }
        return StopReason.MAX_ROUNDS;
    }

    /**
     * The NJ tree of the alignment, rooted and dated (§7.1), in place of the tree as it is.
     *
     * @param first the search's first start, which the log's distances are to; otherwise the second start
     * @return whether the tree is now the NJ tree
     */
    private boolean njStart(boolean first) {
        final String fallback = first ? "; starting from the tree as it is" : "; no second start";
        final dr.evolution.alignment.PatternList patterns = TimeTreeInitializer.patterns(tree, objective);
        if (patterns == null) {
            say("NJ start: no BEAGLE tree likelihood on the tree" + fallback);
            return false;
        }
        final TimeTreeInitializer initializer = new TimeTreeInitializer(tree, rates);
        final double before = evaluator.evaluate();
        if (initializer.adopt(patterns, constraints)) {
            if (first) {
                startTree = new FlexibleTree(tree);
            }
            say("NJ start (" + (initializer.isDated() ? "rooted by root-to-tip regression" : "minimum-variance root") +
                    ", rate " + String.format("%.4g", initializer.getRate()) + "): F = " + format(before) + " -> " +
                    format(evaluator.evaluate()));
            return true;
        }
        say("NJ start: the NJ tree breaks a topology constraint" + fallback);
        return false;
    }

    // the current tree into the candidate set
    private void addCandidate() {
        final double f = evaluator.evaluate();
        if (f > fBest) {
            fBest = f;
            bestState = snapshot();
        }
        candidates.add(TopologyHash.of(tree), f, new FlexibleTree(tree), snapshot(), climb);
    }

    private void finish(TreeStateSnapshot best) {
        if (evaluator.evaluate() < fBest) {
            best.restore();
        }
        // the final fit, with the final maximizer's own settings
        fFinal = refit(finalMaximizer, null, "t3b_final", rounds, evaluator.evaluate());
        fBest = Math.max(fBest, fFinal);
        bestTree = new FlexibleTree(tree);
        nniDelta = settings.nniDelta ? nniDelta() : null;
        seconds = seconds();
        say("stopped (" + stopReason + ") after " + rounds + " rounds: F = " + format(fFinal) + " (" +
                format(fFinal - fStart) + " over the start), " + accepted + " moves accepted, " +
                String.format("%.1fs", seconds));
        try {
            if (settings.logFileName != null) {
                log.writeLog(new File(settings.logFileName));
            }
            if (settings.treeFileName != null) {
                log.writeTree(new File(settings.treeFileName), fFinal, stopReason.toString(), nniDelta);
            }
            if (settings.candidatesFileName != null) {
                log.writeCandidates(new File(settings.candidatesFileName), candidates, bestTree);
            }
        } catch (IOException e) {
            throw new RuntimeException("Could not write the search output: " + e.getMessage(), e);
        }
    }

    /**
     * Local support at the end (§12): at each internal node but the root, F minus the best of its two NNI alternatives,
     * each re-fitted locally (T1), the better also with all heights and rates (T2). An aLRT-like gap; negative when
     * an alternative scores higher, which the search should have taken (counted, and reported).
     *
     * @return by node number, NaN where no alternative is allowed
     */
    private double[] nniDelta() {
        final long started = System.nanoTime();
        final double f = evaluator.evaluate();
        final double[] delta = new double[tree.getNodeCount()];
        java.util.Arrays.fill(delta, Double.NaN);
        nniDeltaNegative = 0;
        for (int i = 0; i < tree.getInternalNodeCount(); ++i) {
            final NodeRef node = tree.getInternalNode(i);
            if (tree.isRoot(node)) {
                continue;
            }
            final BitSet only = new BitSet(tree.getNodeCount());
            only.set(node.getNumber());
            TrialEvaluator.Trial better = null;
            for (RootedNniMove move : NniNeighbourhood.enumerate(tree, constraints, only)) {
                TrialEvaluator.Trial trial = evaluator.trial(move, Double.NEGATIVE_INFINITY, false);
                if (better == null || trial.best > better.best) {
                    better = trial;
                }
            }
            if (better == null) {
                continue;
            }
            final TrialEvaluator.Trial refitted = evaluator.trial(better.move, Double.NEGATIVE_INFINITY, true);
            delta[node.getNumber()] = f - Math.max(better.best, refitted.best);
            if (delta[node.getNumber()] < -settings.moveEpsilon) {
                nniDeltaNegative++;
            }
        }
        trialSeconds += (System.nanoTime() - started) / 1E9;
        if (nniDeltaNegative > 0) {
            say("nniDelta: " + nniDeltaNegative + " node(s) have an NNI alternative scoring higher than the best tree");
        }
        return delta;
    }

    /**
     * One round (§7.3): score the allowed NNIs, near the last change only when speedNNI is on (and all of them before
     * giving up); apply every non-conflicting gaining move at once, falling back to the best single one if the batch
     * ends lower; when none gains, give the closest candidates a T2 look, and the near-ties a longer one.
     *
     * @return the number of moves applied
     */
    private int nniRound(int round) {
        final long started = System.nanoTime();
        final double fInc = evaluator.evaluate();
        final TrialEvaluator.Theta base = evaluator.capture();
        final long t1Before = evaluator.getT1Count();
        final long t2Before = evaluator.getT2Count();

        final BitSet window = settings.speedNni && lastApplied != null
                ? NniNeighbourhood.within(tree, lastApplied, settings.speedNniRadius) : null;
        List<TrialEvaluator.Trial> trials = scan(NniNeighbourhood.enumerate(tree, constraints, window), fInc);
        String scope = window != null ? "local" : "full";
        if (window != null && gaining(trials, fInc).isEmpty() && checkBudgets() == null) {
            // nothing gains near the last change: the rest of the tree before stopping
            Set<RootedNniMove> scored = new HashSet<RootedNniMove>();
            for (TrialEvaluator.Trial trial : trials) {
                scored.add((RootedNniMove) trial.move);
            }
            List<RootedNniMove> rest = new ArrayList<RootedNniMove>();
            for (RootedNniMove move : NniNeighbourhood.enumerate(tree, constraints)) {
                if (!scored.contains(move)) {
                    rest.add(move);
                }
            }
            trials.addAll(scan(rest, fInc));
            scope = "local+full";
        }
        final int nT0 = trials.size();

        final List<TrialEvaluator.Trial> positive = gaining(trials, fInc);
        List<TrialEvaluator.Trial> applied = new ArrayList<TrialEvaluator.Trial>();
        String note = scope;
        double f = fInc;
        if (positive.isEmpty()) {
            TrialEvaluator.Trial winner = checkBudgets() == null ? rescue(trials, fInc) : null;
            if (winner != null) {
                f = evaluator.accept(winner);
                applied.add(winner);
                note += winner.f2 == winner.best ? ",t2" : "";
            }
        } else {
            final TrialEvaluator.Trial single = positive.get(0);
            final List<TrialEvaluator.Trial> batch = settings.batch ? nonConflicting(positive) : null;
            if (batch == null || batch.size() == 1) {
                f = evaluator.accept(single);
                applied.add(single);
            } else {
                final TreeStateSnapshot before = snapshot();
                f = evaluator.acceptBatch(batch, base);
                if (!Double.isNaN(f)) {
                    f = refit(roundMaximizer, evaluator.getT2Budget(), "t2_batch", round, f);
                }
                if (Double.isNaN(f) || f < single.best - settings.batchTolerance) {
                    before.restore();
                    f = evaluator.accept(single);
                    applied.add(single);
                    note += ",batch " + batch.size() + " fell back";
                } else {
                    applied.addAll(batch);
                    note += ",batch " + batch.size();
                }
            }
        }

        if (applied.isEmpty()) {
            lastApplied = null;
        } else {
            List<Integer> touched = new ArrayList<Integer>();
            for (TrialEvaluator.Trial trial : applied) {
                for (int node : trial.move.touchedNodes()) {
                    touched.add(node);
                }
            }
            lastApplied = new int[touched.size()];
            for (int k = 0; k < lastApplied.length; ++k) {
                lastApplied[k] = touched.get(k);
            }
        }

        final int nT1 = (int) (evaluator.getT1Count() - t1Before);
        final int nT2 = (int) (evaluator.getT2Count() - t2Before);
        trialSeconds += (System.nanoTime() - started) / 1E9;
        final TrialEvaluator.Trial first = applied.isEmpty() ? null : applied.get(0);
        log.row(round, "nni", nT0, nT1, nT2, applied.size(), fInc, f, Math.max(fBest, f), describe(applied),
                first != null ? first.f0 - fInc : Double.NaN,
                first != null ? first.f1 - fInc : Double.NaN,
                first != null ? first.f2 - fInc : Double.NaN,
                startTree, seconds(), evaluations(), note, guardReverts);
        say("round " + round + " NNI: F = " + format(f) + " (" + format(f - fInc) + ") " +
                (applied.isEmpty() ? "no move gains" : applied.size() + " move(s): " + describe(applied)) + ", " +
                nT0 + " moves (" + note + "), " + nT1 + " T1, " + nT2 + " T2, " + String.format("%.1fs", seconds()));
        return applied.size();
    }

    private List<TrialEvaluator.Trial> scan(List<? extends TimeTreeMove> moves, double fInc) {
        List<TrialEvaluator.Trial> trials = new ArrayList<TrialEvaluator.Trial>();
        for (TimeTreeMove move : moves) {
            if (checkBudgets() != null) {
                break;
            }
            trials.add(evaluator.trial(move, fInc - settings.deltaT1, false));
        }
        return trials;
    }

    // the trials that gain more than moveEpsilon, best first (ties keep the enumeration order)
    private List<TrialEvaluator.Trial> gaining(List<TrialEvaluator.Trial> trials, double fInc) {
        List<TrialEvaluator.Trial> positive = new ArrayList<TrialEvaluator.Trial>();
        for (TrialEvaluator.Trial trial : trials) {
            if (trial.best > fInc + settings.moveEpsilon) {
                positive.add(trial);
            }
        }
        sortByValue(positive);
        return positive;
    }

    private static void sortByValue(List<TrialEvaluator.Trial> trials) {
        Collections.sort(trials, new Comparator<TrialEvaluator.Trial>() {
            public int compare(TrialEvaluator.Trial x, TrialEvaluator.Trial y) {
                return Double.compare(y.best, x.best);
            }
        });
    }

    /**
     * A greedy batch, best first: a move joins when no node of its quartet is within one edge of a quartet already in
     * the batch. Their local blocks are then disjoint, so their changes can be applied together.
     */
    private List<TrialEvaluator.Trial> nonConflicting(List<TrialEvaluator.Trial> positive) {
        List<TrialEvaluator.Trial> batch = new ArrayList<TrialEvaluator.Trial>();
        BitSet taken = new BitSet(tree.getNodeCount());
        for (TrialEvaluator.Trial trial : positive) {
            final int[] quartet = ((RootedNniMove) trial.move).quartet(tree);
            boolean conflict = false;
            for (int node : quartet) {
                conflict |= taken.get(node);
            }
            if (!conflict) {
                batch.add(trial);
                taken.or(NniNeighbourhood.within(tree, quartet, 1));
            }
        }
        return batch;
    }

    /**
     * The T2 rescue when no move gains: the closest toplistSize candidates within deltaKeep are kept, the best
     * t2Candidates of them re-scored with all heights and rates re-fitted, and those that end within
     * secondChanceWindow below the incumbent re-scored again with a longer T2.
     *
     * @return the best candidate that gains, or null
     */
    private TrialEvaluator.Trial rescue(List<TrialEvaluator.Trial> trials, double fInc) {
        List<TrialEvaluator.Trial> toplist = new ArrayList<TrialEvaluator.Trial>();
        for (TrialEvaluator.Trial trial : trials) {
            if (trial.best > fInc - settings.deltaKeep) {
                toplist.add(trial);
            }
        }
        sortByValue(toplist);
        if (toplist.size() > settings.toplistSize) {
            toplist = new ArrayList<TrialEvaluator.Trial>(toplist.subList(0, settings.toplistSize));
        }

        TrialEvaluator.Trial winner = null;
        List<TrialEvaluator.Trial> nearTies = new ArrayList<TrialEvaluator.Trial>();
        for (int k = 0; k < Math.min(settings.t2Candidates, toplist.size()) && checkBudgets() == null; ++k) {
            TrialEvaluator.Trial trial = evaluator.trial(toplist.get(k).move, Double.NEGATIVE_INFINITY, true);
            if (trial.best > fInc + settings.moveEpsilon) {
                if (winner == null || trial.best > winner.best) {
                    winner = trial;
                }
            } else if (trial.best > fInc - settings.secondChanceWindow) {
                nearTies.add(trial);
            }
        }
        if (winner == null && settings.secondChanceIterations > 0) {
            final MaximizerWrtParameter.Budget longer =
                    new MaximizerWrtParameter.Budget(settings.secondChanceIterations).setQuiet(true);
            for (TrialEvaluator.Trial nearTie : nearTies) {
                if (checkBudgets() != null) {
                    break;
                }
                secondChances++;
                TrialEvaluator.Trial trial = evaluator.trial(nearTie.move, Double.NEGATIVE_INFINITY, longer);
                if (trial.best > fInc + settings.moveEpsilon && (winner == null || trial.best > winner.best)) {
                    winner = trial;
                }
            }
        }
        return winner;
    }

    /**
     * The root scan of §5.2, after a new best: every root position one or two edges away (root-child NNIs, alone and
     * chained), each with T0 and T1 on its block, and T2 on the best of them if none gains yet.
     *
     * @return whether the root moved
     */
    private boolean rootScan(int round) {
        final long started = System.nanoTime();
        final double fInc = evaluator.evaluate();
        final long t1Before = evaluator.getT1Count();
        final long t2Before = evaluator.getT2Count();
        final List<TrialEvaluator.Trial> trials = scan(NniNeighbourhood.rootMoves(tree, constraints), fInc);
        rootScans++;

        TrialEvaluator.Trial winner = null;
        for (TrialEvaluator.Trial trial : trials) {
            if (winner == null || trial.best > winner.best) {
                winner = trial;
            }
        }
        if (winner != null && !(winner.best > fInc + settings.moveEpsilon) && checkBudgets() == null) {
            winner = evaluator.trial(winner.move, Double.NEGATIVE_INFINITY, true);
        }
        if (winner != null && !(winner.best > fInc + settings.moveEpsilon)) {
            winner = null;
        }
        double f = fInc;
        if (winner != null) {
            f = evaluator.accept(winner);
            lastApplied = winner.move.touchedNodes();
            rootMoves++;
        }
        trialSeconds += (System.nanoTime() - started) / 1E9;
        log.row(round, "root", trials.size(), (int) (evaluator.getT1Count() - t1Before),
                (int) (evaluator.getT2Count() - t2Before), winner != null ? 1 : 0, fInc, f, Math.max(fBest, f),
                winner != null ? winner.move.describe() : null,
                winner != null ? winner.f0 - fInc : Double.NaN,
                winner != null ? winner.f1 - fInc : Double.NaN,
                winner != null ? winner.f2 - fInc : Double.NaN,
                startTree, seconds(), evaluations(), "", guardReverts);
        say("round " + round + " root scan: F = " + format(f) + " (" + format(f - fInc) + ") " +
                (winner != null ? winner.move.describe() : "the root stays") + ", " + trials.size() + " positions");
        return winner != null;
    }

    /**
     * The SPR phase (plan 7 §7.3, after RAxML-NG): fast rounds within a radius found by trying 5, 10, ... until the
     * gain stops, repeated while they gain more than 1; then thorough rounds over windows of radii [1, 5], [6, 10],
     * ..., back to [1, 5] after a gain.
     *
     * @return the number of moves accepted
     */
    private int sprPhase(int round) {
        final int maxRadius = Math.min(settings.sprMaxRadius, tree.getExternalNodeCount() - 3);
        if (maxRadius < 1) {
            return 0;
        }
        final int acceptedBefore = sprAccepted;
        if (sprCutoff < 0) {
            sprCutoff = Math.abs(evaluator.evaluate()) / 1000.0;
        }
        int radius = settings.sprRadius > 0 ? Math.min(settings.sprRadius, maxRadius) : 0;
        if (radius == 0) {
            // auto-detect: widen while widening still gains
            radius = Math.min(5, maxRadius);
            while (checkBudgets() == null) {
                double gain = sprRound(round, 1, radius, false);
                if (gain < settings.epsilon || radius >= maxRadius) {
                    break;
                }
                radius = Math.min(radius + 5, maxRadius);
            }
        }
        while (checkBudgets() == null && sprRound(round, 1, radius, false) > 1.0) {
            // fast rounds while they gain
        }
        int rMin = 1;
        while (rMin <= maxRadius && checkBudgets() == null) {
            final int rMax = Math.min(rMin + 4, maxRadius);
            final double gain = sprRound(round, rMin, rMax, true);
            rMin = gain > settings.epsilon ? 1 : rMax + 1;
        }
        return sprAccepted - acceptedBefore;
    }

    /**
     * One SPR round over the radii [rMin, rMax]. Each prunable subtree, in priority order, gets its targets scored
     * breadth-first at T0, a target's neighbours visited only while its score is within the cutoff of the incumbent;
     * in thorough rounds the best few also get T1. The best target is accepted at once if it gains (greedy, first
     * improvement), else kept on a toplist if within deltaKeep. At the end the toplist's moves still valid are
     * re-scored at T1, the best t2Candidates of them at T2, and the best accepted if it gains; the cutoff becomes
     * subtreeCutoff times the mean drop seen, and a round with moves ends with T3a.
     *
     * @return F after the round minus F before it
     */
    private double sprRound(int round, int rMin, int rMax, boolean thorough) {
        final long started = System.nanoTime();
        final double fStart = evaluator.evaluate();
        final long t1Before = evaluator.getT1Count();
        final long t2Before = evaluator.getT2Count();
        final int toplistSize = thorough ? settings.toplistSize : settings.sprFastToplistSize;
        final List<TrialEvaluator.Trial> toplist = new ArrayList<TrialEvaluator.Trial>();
        final double[] drops = new double[2]; // sum and count of T0 drops below the incumbent
        final int[] scored = new int[1];
        int moves = 0;
        double fInc = fStart;
        sprRounds++;

        for (int i : SprNeighbourhood.pruneOrder(tree)) {
            if (checkBudgets() != null) {
                break;
            }
            final NodeRef node = tree.getNode(i);
            if (tree.isRoot(node) || tree.isRoot(tree.getParent(node))) {
                continue; // an accepted move has moved it next to the root
            }
            final double incumbent = fInc;
            final List<TrialEvaluator.Trial> trials = new ArrayList<TrialEvaluator.Trial>();
            SprNeighbourhood.search(tree, i, rMin, rMax, constraints, new SprNeighbourhood.Visitor() {
                public double visit(ClampedSprMove move, int distance) {
                    TrialEvaluator.Trial trial = evaluator.trial(move, Double.POSITIVE_INFINITY, false);
                    trials.add(trial);
                    scored[0]++;
                    if (trial.f0 < incumbent && trial.f0 > Double.NEGATIVE_INFINITY) {
                        drops[0] += incumbent - trial.f0;
                        drops[1]++;
                    }
                    return trial.f0;
                }
            }, incumbent - sprCutoff);
            if (trials.isEmpty()) {
                continue;
            }
            sortByValue(trials);
            TrialEvaluator.Trial best = trials.get(0);
            if (thorough) {
                // T1 on the best target and up to three more within deltaT1 of it
                for (int k = 0; k < Math.min(4, trials.size()) && trials.get(k).f0 > best.f0 - settings.deltaT1; ++k) {
                    TrialEvaluator.Trial refined = evaluator.trial(trials.get(k).move, Double.NEGATIVE_INFINITY, false);
                    if (k == 0 || refined.best > best.best) {
                        best = refined;
                    }
                }
            }
            if (best.best > fInc + settings.moveEpsilon) {
                fInc = evaluator.accept(best);
                moves++;
            } else if (best.best > fInc - settings.deltaKeep) {
                toplist.add(best);
                sortByValue(toplist);
                if (toplist.size() > toplistSize) {
                    toplist.remove(toplist.size() - 1);
                }
            }
        }

        // the toplist: T1 on the moves still valid, T2 on the best of them
        if (checkBudgets() == null) {
            final List<TrialEvaluator.Trial> refined = new ArrayList<TrialEvaluator.Trial>();
            for (TrialEvaluator.Trial kept : toplist) {
                if (kept.move.isValid(tree) && ((ClampedSprMove) kept.move).isFeasible(tree)
                        && !violates((ClampedSprMove) kept.move) && checkBudgets() == null) {
                    refined.add(evaluator.trial(kept.move, Double.NEGATIVE_INFINITY, false));
                }
            }
            sortByValue(refined);
            TrialEvaluator.Trial winner = null;
            for (int k = 0; k < Math.min(settings.t2Candidates, refined.size()) && checkBudgets() == null; ++k) {
                TrialEvaluator.Trial trial = evaluator.trial(refined.get(k).move, Double.NEGATIVE_INFINITY, true);
                if (trial.best > fInc + settings.moveEpsilon && (winner == null || trial.best > winner.best)) {
                    winner = trial;
                }
            }
            if (winner != null) {
                fInc = evaluator.accept(winner);
                moves++;
            }
        }

        if (drops[1] > 0) {
            sprCutoff = Math.max(1.0, settings.subtreeCutoff * drops[0] / drops[1]);
        }
        sprAccepted += moves;
        trialSeconds += (System.nanoTime() - started) / 1E9;
        final String phase = thorough ? "spr_slow" : "spr_fast";
        log.row(round, phase, scored[0], (int) (evaluator.getT1Count() - t1Before),
                (int) (evaluator.getT2Count() - t2Before), moves, fStart, fInc, Math.max(fBest, fInc), null,
                Double.NaN, Double.NaN, Double.NaN, startTree, seconds(), evaluations(),
                "r " + rMin + "-" + rMax + ", cutoff " + String.format("%.2f", sprCutoff), guardReverts);
        say("round " + round + " SPR " + (thorough ? "slow" : "fast") + " r=" + rMin + "-" + rMax + ": F = " +
                format(fInc) + " (" + format(fInc - fStart) + "), " + moves + " moves, " + scored[0] + " targets, " +
                String.format("%.1fs", seconds()));
        double f = fInc;
        if (moves > 0) {
            f = refit(roundMaximizer, t3aBudget(), "t3a", round, evaluator.evaluate());
        }
        return f - fStart;
    }

    /**
     * One round of subtree swaps (BEAST's wide exchange) whose parents are 2 to swapMaxDistance nodes apart.
     *
     * @return the number of moves applied (0 or 1)
     */
    private int swapRound(int round) {
        swapRounds++;
        final int applied = oneMoveRound(round, "swap", SubtreeSwapMove.enumerate(tree, 2, settings.swapMaxDistance,
                constraints), "d 2-" + settings.swapMaxDistance);
        swapAccepted += applied;
        return applied;
    }

    /**
     * One round of root-changing SPR (§5.4): every subtree above the root, and each root child with the root to the
     * branches at most rootMaxDepth below its sibling. The root gets another node; node numbers stay with their nodes,
     * and snapshots taken before restore by clade.
     *
     * @return the number of moves applied (0 or 1)
     */
    private int rootSprRound(int round) {
        rootSprRounds++;
        final int applied = oneMoveRound(round, "root_spr", RootChangingSprMove.enumerate(tree, settings.rootMaxDepth,
                constraints), "depth " + settings.rootMaxDepth);
        rootSprAccepted += applied;
        return applied;
    }

    /**
     * T0 and T1 on every move, the best accepted if it gains, else T2 on the closest t2Candidates within deltaKeep.
     *
     * @return the number of moves applied (0 or 1)
     */
    private int oneMoveRound(int round, String phase, List<? extends TimeTreeMove> moves, String note) {
        final long started = System.nanoTime();
        final double fInc = evaluator.evaluate();
        final long t1Before = evaluator.getT1Count();
        final long t2Before = evaluator.getT2Count();
        final List<TrialEvaluator.Trial> trials = scan(moves, fInc);
        final List<TrialEvaluator.Trial> positive = gaining(trials, fInc);
        TrialEvaluator.Trial winner = positive.isEmpty() ? null : positive.get(0);
        if (winner == null && checkBudgets() == null) {
            List<TrialEvaluator.Trial> close = new ArrayList<TrialEvaluator.Trial>();
            for (TrialEvaluator.Trial trial : trials) {
                if (trial.best > fInc - settings.deltaKeep) {
                    close.add(trial);
                }
            }
            sortByValue(close);
            for (int k = 0; k < Math.min(settings.t2Candidates, close.size()) && checkBudgets() == null; ++k) {
                TrialEvaluator.Trial trial = evaluator.trial(close.get(k).move, Double.NEGATIVE_INFINITY, true);
                if (trial.best > fInc + settings.moveEpsilon && (winner == null || trial.best > winner.best)) {
                    winner = trial;
                }
            }
        }
        double f = fInc;
        if (winner != null) {
            f = evaluator.accept(winner);
        }
        trialSeconds += (System.nanoTime() - started) / 1E9;
        log.row(round, phase, trials.size(), (int) (evaluator.getT1Count() - t1Before),
                (int) (evaluator.getT2Count() - t2Before), winner != null ? 1 : 0, fInc, f, Math.max(fBest, f),
                winner != null ? winner.move.describe() : null,
                winner != null ? winner.f0 - fInc : Double.NaN,
                winner != null ? winner.f1 - fInc : Double.NaN,
                winner != null ? winner.f2 - fInc : Double.NaN,
                startTree, seconds(), evaluations(), note, guardReverts);
        say("round " + round + " " + phase + ": F = " + format(f) + " (" + format(f - fInc) + ") " +
                (winner != null ? winner.move.describe() : "no move gains") + ", " + trials.size() + " moves");
        return winner != null ? 1 : 0;
    }

    private boolean violates(ClampedSprMove move) {
        return constraints.size() > 0 && constraints.violates(tree, move, TopologyConstraints.clades(tree));
    }

    private static String describe(List<TrialEvaluator.Trial> applied) {
        if (applied.isEmpty()) {
            return null;
        }
        StringBuilder sb = new StringBuilder();
        for (TrialEvaluator.Trial trial : applied) {
            sb.append(sb.length() > 0 ? "; " : "").append(trial.move.describe());
        }
        return sb.toString();
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
        } else if (phase.startsWith("t3b")) {
            t3bSeconds += elapsed;
        } // a batch's T2 is inside its round's trial time
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

    public int getSecondChances() {
        return secondChances;
    }

    public int getPerturbations() {
        return perturbations;
    }

    public CandidateSet getCandidates() {
        return candidates;
    }

    /**
     * @return by node number, the final F minus the best NNI alternative at the node (NaN for none), or null
     */
    public double[] getNniDelta() {
        return nniDelta;
    }

    public int getSprRounds() {
        return sprRounds;
    }

    public int getSwapRounds() {
        return swapRounds;
    }

    public int getSwapAcceptedMoves() {
        return swapAccepted;
    }

    public int getRootSprRounds() {
        return rootSprRounds;
    }

    public int getRootSprAcceptedMoves() {
        return rootSprAccepted;
    }

    public int getSprAcceptedMoves() {
        return sprAccepted;
    }

    public int getRootScans() {
        return rootScans;
    }

    public int getRootMoves() {
        return rootMoves;
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

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.FRAMEWORK;
    }

    @Override
    public String getDescription() {
        return "Starting the topology search from the neighbour-joining tree";
    }

    // cited only when a start is the NJ tree
    @Override
    public List<Citation> getCitations() {
        return settings.njStart ? Collections.singletonList(TimeTreeInitializer.NJ_CITATION) :
                Collections.<Citation>emptyList();
    }

    public SearchSettings getSettings() {
        return settings;
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
        StringBuilder sb = new StringBuilder("Topology search (rooted NNI" + (settings.spr ? ", SPR" : "") +
                (settings.swap ? ", subtree swap" : "") + " hill-climb)\n");
        if (!executed) {
            return sb.append("Not yet executed.\n").toString();
        }
        sb.append("Settings: ").append(settings.describe()).append("\n");
        sb.append("StopReason: ").append(stopReason).append("\n");
        sb.append("Rounds: ").append(rounds).append("\n");
        sb.append("AcceptedMoves: ").append(accepted).append("\n");
        sb.append("Trials: ").append(evaluator.getTrials()).append(" (T1 ").append(evaluator.getT1Count())
                .append(", T2 ").append(evaluator.getT2Count()).append(")\n");
        sb.append("SecondChances: ").append(secondChances).append("\n");
        if (perturbation != null) {
            sb.append("Perturbations: ").append(perturbations).append(" (operators ")
                    .append(perturbation.getOperatorNames()).append(", ").append(perturbation.getKept())
                    .append(" of ").append(perturbation.getDraws()).append(" draws kept)\n");
            sb.append("Candidates: ").append(candidates.size()).append(", best from climb ")
                    .append(candidates.size() > 0 ? candidates.getCandidates().get(0).climb : 0).append("\n");
        }
        if (nniDelta != null) {
            double min = Double.POSITIVE_INFINITY;
            for (double d : nniDelta) {
                if (!Double.isNaN(d)) {
                    min = Math.min(min, d);
                }
            }
            sb.append("NniDelta: smallest ").append(min).append(nniDeltaNegative > 0 ? " (" + nniDeltaNegative +
                    " negative: an NNI scores higher)" : "").append("\n");
        }
        if (settings.swap) {
            sb.append("SwapRounds: ").append(swapRounds).append(" (moves accepted ").append(swapAccepted).append(")\n");
        }
        if (settings.rootMoves) {
            sb.append("RootSprRounds: ").append(rootSprRounds).append(" (moves accepted ").append(rootSprAccepted)
                    .append(")\n");
        }
        if (settings.spr) {
            sb.append("SprRounds: ").append(sprRounds).append(" (moves accepted ").append(sprAccepted).append(")\n");
        }
        sb.append("RootScans: ").append(rootScans).append(" (root moved ").append(rootMoves).append(")\n");
        sb.append("GuardReverts: ").append(guardReverts).append("\n");
        if (evaluator.getUnderflowRecomputes() > 0) {
            sb.append("UnderflowRecomputes: ").append(evaluator.getUnderflowRecomputes())
                    .append(" (BEAGLE rescaled inside a trial; everything recomputed after it)\n");
        }
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
