/*
 * ExplicitHitsHistoryModel.java
 *
 * Copyright © 2026, the BEAST Development Team.
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
 */

package dr.evomodel.treedatalikelihood.MAT;

import dr.math.MathUtils;

/**
 * Explicit hits: each site's history on a branch is a path of states a = x_0 -> x_1 -> ... -> x_k = b
 * (any k, including back-mutations when a == b), stored as the site's substitutions in order. The
 * times of the hits are integrated out, so a path contributes
 * sum_i log(nu Q_{x_{i-1} x_i}) + log int_{0 < t_1 < ... < t_k < L} exp(-sum_i nu q_{x_i} (t_{i+1} - t_i)).
 * Summed over paths this gives the transition probability P_ab(nu L), so the ancestral states are
 * drawn as for integrated hits and then each site's path is drawn given its end states.
 */
public class ExplicitHitsHistoryModel extends AbstractBranchHistoryModel {

    public static final String NAME = "explicitHits";

    private final int[][] pathStates = new int[1][];
    private final double[][] pathTimes = new double[1][];
    private double[] rates = new double[8];

    public ExplicitHitsHistoryModel(int stateCount) {
        super(stateCount);
    }

    @Override
    public String getName() {
        return NAME;
    }

    @Override
    public boolean hasExplicitHistories() {
        return true;
    }

    @Override
    public void getLogTransitionDensities(double length, double nu, double[] out) {
        System.arraycopy(getLogTransitionProbabilities(nu * length), 0, out, 0, stateCount * stateCount);
    }

    @Override
    public double getBranchLogLikelihood(BranchMutations mutations, double length, double exitRateSum, double[] composition,
                                         double[] patternWeights, double[] patternRates) {
        double term = -exitRateSum * length;
        int start = 0;
        while (start < mutations.getSubstitutionCount()) {
            int end = mutations.getSiteEnd(start);
            int p = mutations.getPattern(start);
            double nu = patternRates[p];
            // the -lambda L term has already counted the site in its final state over the whole branch
            term += patternWeights[p] * (nu * exitRates[mutations.getToState(end - 1)] * length
                    + getSiteLogLikelihood(mutations, start, end, length, nu));
            start = end;
        }
        return term;
    }

    @Override
    public double getUnchangedSiteLogLikelihood(int a, double length, double nu) {
        return -nu * exitRates[a] * length;
    }

    @Override
    public double getSiteLogLikelihood(BranchMutations mutations, int start, int end, double length, double nu) {
        int k = end - start;
        int[] states = new int[k + 1];
        double[] times = new double[k + 1];
        states[0] = mutations.getFromState(start);
        for (int i = 0; i < k; i++) {
            states[i + 1] = mutations.getToState(start + i);
            times[i + 1] = mutations.hasTimes() ? mutations.getTime(start + i) : (i + 1.0) / (k + 1.0);
        }
        return pathLogDensity(states, times, k, length, nu);
    }

    /**
     * @return the log density of a path of k hits (times as fractions of the branch, times[1..k])
     */
    protected double pathLogDensity(int[] states, double[] times, int k, double length, double nu) {
        if (rates.length < k + 1) {
            rates = new double[2 * (k + 1)];
        }
        double logRate = 0.0;
        for (int i = 0; i <= k; i++) {
            rates[i] = nu * exitRates[states[i]];
            if (i > 0) {
                logRate += Math.log(nu * infinitesimalMatrix[states[i - 1] * stateCount + states[i]]);
            }
        }
        return logRate + logOrderedTimesIntegral(rates, k, length);
    }

    @Override
    public double getHistoryLogProposal(BranchMutations mutations, int start, int end, int fromState, int toState,
                                        double length, double nu) {
        double logDensity = start == end ? getUnchangedSiteLogLikelihood(fromState, length, nu) :
                getSiteLogLikelihood(mutations, start, end, length, nu);
        return logDensity - getLogTransitionProbabilities(nu * length)[fromState * stateCount + toState];
    }

    @Override
    public double sampleHistory(int pattern, int fromState, int toState, double length, double nu, BranchMutations.Builder builder) {
        double logP = getLogTransitionProbabilities(nu * length)[fromState * stateCount + toState];
        if (fromState == toState) {
            // the empty path, or (conditional on not being empty) a path of back-mutations
            double logEmpty = getUnchangedSiteLogLikelihood(fromState, length, nu) - logP;
            if (MathUtils.nextDouble() < Math.exp(logEmpty)) {
                return logEmpty;
            }
        }
        int k;
        do {
            k = samplePath(fromState, toState, length, nu, pathStates, pathTimes);
        } while (k == 0);
        for (int i = 1; i <= k; i++) {
            builder.addSubstitution(pattern, pathStates[0][i - 1], pathStates[0][i], pathTimes[0][i]);
        }
        return pathLogDensity(pathStates[0], pathTimes[0], k, length, nu) - logP;
    }
}
