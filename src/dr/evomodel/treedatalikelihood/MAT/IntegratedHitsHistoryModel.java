/*
 * IntegratedHitsHistoryModel.java
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
 * Any number of hits per site per branch, summed analytically: a site going from state a to state b
 * along a branch contributes log P_ab(nu L), the transition probability of the substitution model.
 * Summed over the ancestral states this is the standard (Felsenstein) likelihood. The hits themselves
 * are not represented; the branch's events are the differences between its end states.
 *
 * The branch term is sum_a n_a log P_aa(L) plus a correction for each substitution, where n_a is the
 * summed weight of the node's sites in state a, so it needs the same rate at every site (no per-site
 * rates).
 */
public class IntegratedHitsHistoryModel extends AbstractBranchHistoryModel {

    public static final String NAME = "integratedHits";

    public IntegratedHitsHistoryModel(int stateCount) {
        super(stateCount);
    }

    @Override
    public String getName() {
        return NAME;
    }

    @Override
    public void checkCompatibility(MATDataLikelihoodDelegate delegate) {
        if (delegate.hasSiteRates()) {
            throw new IllegalArgumentException("The " + NAME + " branch history model does not support per-site rates");
        }
    }

    @Override
    public boolean hasExplicitHistories() {
        return false;
    }

    @Override
    public boolean needsComposition() {
        return true;
    }

    @Override
    public boolean needsExitRate() {
        return false;
    }

    @Override
    public void getLogTransitionDensities(double length, double nu, double[] out) {
        System.arraycopy(getLogTransitionProbabilities(nu * length), 0, out, 0, stateCount * stateCount);
    }

    @Override
    public double getBranchLogLikelihood(BranchMutations mutations, double length, double exitRateSum, double[] composition,
                                         double[] patternWeights, double[] patternRates) {
        // every site has the same rate
        double nu = patternRates.length > 0 ? patternRates[0] : 1.0;
        double[] logP = getLogTransitionProbabilities(nu * length);
        double term = 0.0;
        for (int a = 0; a < stateCount; a++) {
            if (composition[a] > 0.0) {
                term += composition[a] * logP[a * stateCount + a];
            }
        }
        for (int k = 0; k < mutations.getSubstitutionCount(); k++) {
            int p = mutations.getPattern(k);
            int a = mutations.getFromState(k);
            int b = mutations.getToState(k);
            // the composition counted the site as staying in state b
            term += patternWeights[p] * (logP[a * stateCount + b] - logP[b * stateCount + b]);
        }
        return term;
    }

    @Override
    public double getUnchangedSiteLogLikelihood(int a, double length, double nu) {
        return getLogTransitionProbabilities(nu * length)[a * stateCount + a];
    }

    @Override
    public double getSiteLogLikelihood(BranchMutations mutations, int start, int end, double length, double nu) {
        return getLogTransitionProbabilities(nu * length)[mutations.getFromState(start) * stateCount + mutations.getToState(end - 1)];
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
