/*
 * SingleHitHistoryModel.java
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

/**
 * At most one substitution per site per branch, with its time integrated out: a site staying in
 * state a contributes -nu q_a L, and a substitution a->b contributes
 * log(nu Q_ab L) - nu q_b L + log(expm1(x)/x), x = nu (q_b - q_a) L.
 */
public class SingleHitHistoryModel extends AbstractBranchHistoryModel {

    public static final String NAME = "singleHit";

    public SingleHitHistoryModel(int stateCount) {
        super(stateCount);
    }

    @Override
    public String getName() {
        return NAME;
    }

    @Override
    public boolean hasExplicitHistories() {
        return false;
    }

    @Override
    public void getLogTransitionDensities(double length, double nu, double[] out) {
        for (int a = 0; a < stateCount; a++) {
            for (int b = 0; b < stateCount; b++) {
                out[a * stateCount + b] = a == b ? -nu * exitRates[a] * length : substitutionLogDensity(a, b, length, nu);
            }
        }
    }

    private double substitutionLogDensity(int a, int b, double length, double nu) {
        return Math.log(nu * infinitesimalMatrix[a * stateCount + b] * length) - nu * exitRates[b] * length
                + MATDataLikelihoodDelegate.logExpm1OverX(nu * (exitRates[b] - exitRates[a]) * length);
    }

    @Override
    public double getBranchLogLikelihood(BranchMutations mutations, double length, double exitRateSum, double[] composition,
                                         double[] patternWeights, double[] patternRates) {
        double term = -exitRateSum * length;
        for (int k = 0; k < mutations.getSubstitutionCount(); k++) {
            int p = mutations.getPattern(k);
            int a = mutations.getFromState(k);
            int b = mutations.getToState(k);
            double nu = patternRates[p];
            // the -lambda L term has already counted the site in state b over the whole branch
            term += patternWeights[p] * (substitutionLogDensity(a, b, length, nu) + nu * exitRates[b] * length);
        }
        return term;
    }

    @Override
    public double getUnchangedSiteLogLikelihood(int a, double length, double nu) {
        return -nu * exitRates[a] * length;
    }

    @Override
    public double getSiteLogLikelihood(BranchMutations mutations, int start, int end, double length, double nu) {
        return substitutionLogDensity(mutations.getFromState(start), mutations.getToState(start), length, nu);
    }

    @Override
    public double getHistoryLogProposal(BranchMutations mutations, int start, int end, int fromState, int toState,
                                        double length, double nu) {
        return 0.0;
    }

    @Override
    public double sampleHistory(int pattern, int fromState, int toState, double length, double nu, BranchMutations.Builder builder) {
        if (fromState != toState) {
            builder.addSubstitution(pattern, fromState, toState, 0.0);
        }
        return 0.0;
    }
}
