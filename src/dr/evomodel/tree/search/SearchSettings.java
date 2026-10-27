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
    // T1 for candidates whose T0 is within this of the incumbent. On Dengue4 with the HMC relaxed clock every NNI's T0
    // was more than 20 below the incumbent although some gained after T1, and T1 on all of them took 0.3s: so T1 on
    // nearly every move until the tier diagnostics (plan 7 §6.2) tune it
    public double deltaT1 = 1000.0;
    public double deltaKeep = 100.0; // T2 rescue for candidates within this of the incumbent
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

    // budgets (§6.4); 0: no limit
    public int maxRounds = 0; // 0: the number of internal nodes
    public double maxRunTime = 0.0; // seconds
    public long maxEvaluations = 0;

    public BranchRateWarmStart.Mode rateWarmStart = BranchRateWarmStart.Mode.PRESERVE_SUBSTITUTIONS;
    public int verifyEvery = 0;
    public long seed = -1; // < 0: leave the random number generator alone
    public boolean printToScreen = true;
    public String treeFileName = null;
    public String logFileName = null;

    public String describe() {
        return "epsilon = " + epsilon + ", moveEpsilon = " + moveEpsilon + ", deltaT1 = " + deltaT1 +
                ", deltaKeep = " + deltaKeep + ", t2Candidates = " + t2Candidates + ", t1Sweeps = " + t1Sweeps +
                ", rateWarmStart = " + rateWarmStart.getName() + ", maxRounds = " +
                (maxRounds > 0 ? maxRounds : "n - 1") + ", maxRunTime = " + (maxRunTime > 0 ? maxRunTime + "s" : "none") +
                ", maxEvaluations = " + (maxEvaluations > 0 ? maxEvaluations : "none") +
                ", verifyEvery = " + verifyEvery;
    }
}
