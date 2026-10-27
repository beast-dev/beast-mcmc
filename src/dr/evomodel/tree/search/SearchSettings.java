/*
 * SearchSettings.java
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

/**
 * The settings of a topology search (plan 7 §4.3, §6, §7.6), with the Standard preset's values as defaults.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public final class SearchSettings {

    // acceptance (§3.3)
    public double epsilon = 0.1; // a new global best, and a round's gain worth another round
    public double moveEpsilon = 0.01; // a move's gain over the incumbent to accept it

    // tiers (§6.1)
    // T1 for candidates whose T0 is within this of the incumbent, and the T2 rescue for those whose T1 is within
    // deltaKeep. The tier diagnostics (plan 7 §6.2, TierFidelityDiagnostics) on Dengue4: the moves that gain after T2
    // were at most 30.4 below the incumbent at T0 and 29.0 at T1 (strict clock; 4.9 and 3.6 with the HMC clock and
    // rates kept); 50 is the plan's ceiling for an adaptive deltaKeep
    public double deltaT1 = 50.0;
    public double deltaKeep = 50.0;
    public int t1Sweeps = 2;
    public int t2Candidates = 3;
    public int t2Iterations = 20;
    public int t3aIterations = 200;
    public double t3aRelTol = 1E-6;
    public int t3aStallWindow = 10;
    public int t3bStartIterations = 100;
    public double t3bStartRelTol = 1E-5;
    public int t3bBestIterations = 500;
    public double t3bBestRelTol = 1E-7;
    public double modelEpsilon = 0.1; // a new best by more than this re-fits everything (T3b)

    // the NNI tier set (§7.3, M1b)
    public boolean speedNni = true; // after a change, score only the NNIs near it; a full scan before stopping
    public int speedNniRadius = 2; // in edges from the nodes of the moves just applied
    public boolean batch = true; // apply every non-conflicting gaining NNI at once
    public double batchTolerance = 1E-3; // fall back to the best single move if the batch ends lower by this
    public int toplistSize = 20; // the closest candidates of a round, kept for the T2 rescue
    public double secondChanceWindow = 2.0; // a T2 candidate this close below the incumbent gets a longer T2
    public int secondChanceIterations = 80;
    public boolean rootScan = true; // after a new best, the root positions one and two edges away

    // SPR (§5.3, §7.3, M4)
    public boolean spr = false; // when no NNI gains, clamped-height SPR rounds
    public int sprMaxRadius = 10; // the largest distance of a regraft from the pruned subtree's old place
    public int sprRadius = 0; // the fast rounds' radius; 0: found by trying 5, 10, ...
    public double subtreeCutoff = 1.0; // the scan's cutoff is this times the mean drop of the round before
    public int sprFastToplistSize = 60; // the fast rounds' toplist (thorough rounds use toplistSize)

    // subtree swaps (BEAST's wide exchange), after SPR
    public boolean swap = false;
    public int swapMaxDistance = 4; // the farthest apart the two parents may be, in nodes (GibbsSubtreeSwap's 4)

    // root-changing SPR (§5.4, M6), after SPR and swaps: a subtree above the root, or a root child with the root
    public boolean rootMoves = false;
    public int rootMaxDepth = 3; // prune at the root: the deepest target below the new root, in branches

    // the start (§7.1, M6): false starts from the tree as it is; true from the NJ tree, rooted and dated
    public boolean njStart = false;

    // perturbation and the candidate set (§7.2, §7.4, M5)
    public int nStop = 0; // stop after this many perturbed hill-climbs in a row without a new best; 0: no perturbation
    public double perturbStrength = 0.3; // operator moves per perturbation: max(1, floor(strength (n - 2)))
    public String perturbOperators = "nni,wide,leap,wb"; // BEAST operators drawn from, by doOperation
    public double leapSize = 0.1; // SubtreeLeap's distance scale, relative to the root height
    public int candidateSetSize = 20;
    public int parentPool = 5; // perturb a tree drawn uniformly from this many best candidates
    public boolean nniDelta = true; // at the end, F minus the best NNI alternative at each node
    public String candidatesFileName = null;

    // budgets (§6.4); 0: no limit
    public int maxRounds = 0; // 0: the number of internal nodes
    public double maxRunTime = 0.0; // seconds
    public long maxEvaluations = 0;

    // rates stay with their nodes: after a MAP fit the HMC clock's scale is small, so rates rescaled to preserve
    // substitutions cost millions in the rate prior (median T0 -7.5e6 on Dengue4), and T0 and T1 stop ranking moves
    // (Spearman 0.65 against T2, 0.93-0.95 with rates kept; the search ends at -3928.0 instead of -3883.1)
    public BranchRateWarmStart.Mode rateWarmStart = BranchRateWarmStart.Mode.KEEP;
    public int verifyEvery = 0;
    public long seed = -1; // < 0: leave the random number generator alone
    public boolean printToScreen = true;
    public String treeFileName = null;
    public String logFileName = null;

    public String describe() {
        return "epsilon = " + epsilon + ", moveEpsilon = " + moveEpsilon + ", deltaT1 = " + deltaT1 +
                ", deltaKeep = " + deltaKeep + ", t2Candidates = " + t2Candidates + ", t1Sweeps = " + t1Sweeps +
                ", speedNni = " + speedNni + ", batch = " + batch + ", toplistSize = " + toplistSize +
                ", secondChanceWindow = " + secondChanceWindow + ", rootScan = " + rootScan + ", moves = nni" + (spr ? ",spr" : "") + (swap ? ",swap" : "") + (rootMoves ? ",root" : "") +
                (spr ? ", sprMaxRadius = " + sprMaxRadius + ", sprRadius = " + (sprRadius > 0 ? sprRadius : "auto") +
                        ", subtreeCutoff = " + subtreeCutoff : "") +
                ", rateWarmStart = " + rateWarmStart.getName() + ", maxRounds = " +
                (maxRounds > 0 ? maxRounds : "n - 1") + ", maxRunTime = " + (maxRunTime > 0 ? maxRunTime + "s" : "none") +
                ", maxEvaluations = " + (maxEvaluations > 0 ? maxEvaluations : "none") +
                ", verifyEvery = " + verifyEvery + (njStart ? ", startingTree = nj" : "") +
                (nStop > 0 ? ", nStop = " + nStop + ", perturbStrength = " + perturbStrength + ", perturbOperators = " +
                        perturbOperators + ", candidateSetSize = " + candidateSetSize + ", parentPool = " + parentPool : "");
    }
}
