/*
 * TimedHitsHistoryModel.java
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
 * Explicit hits with their positions along the branch (as in Varilly et al., Nature 2026): a path of k hits at fractions
 * u_1 < ... < u_k of the branch contributes
 * sum_i log(nu Q_{x_{i-1} x_i}) - sum_i nu q_{x_i} (u_{i+1} - u_i) L + k log L,
 * the CTMC path density with the change of variables from times to fractions (so that changing the
 * node heights rescales the hits along their branches).
 */
public class TimedHitsHistoryModel extends ExplicitHitsHistoryModel {

    public static final String NAME = "timedHits";

    public TimedHitsHistoryModel(int stateCount) {
        super(stateCount);
    }

    @Override
    public String getName() {
        return NAME;
    }

    @Override
    protected double pathLogDensity(int[] states, double[] times, int k, double length, double nu) {
        double logDensity = 0.0;
        double previous = 0.0;
        for (int i = 0; i <= k; i++) {
            double next = i < k ? times[i + 1] : 1.0;
            logDensity -= nu * exitRates[states[i]] * (next - previous) * length;
            if (i > 0) {
                logDensity += Math.log(nu * infinitesimalMatrix[states[i - 1] * stateCount + states[i]] * length);
            }
            previous = next;
        }
        return logDensity;
    }

    @Override
    public void initialise(DefaultMutationStore store) {
        // hits built by parsimony have no positions: place each site's hits evenly along the branch
        for (int n = 0; n < store.getNodeCount(); n++) {
            BranchMutations mutations = store.getBranchMutations(n);
            if (mutations.getSubstitutionCount() > 0 && !mutations.hasTimes()) {
                BranchMutations.Builder builder = new BranchMutations.Builder(true);
                int start = 0;
                while (start < mutations.getSubstitutionCount()) {
                    int end = mutations.getSiteEnd(start);
                    for (int i = start; i < end; i++) {
                        builder.addSubstitution(mutations.getPattern(i), mutations.getFromState(i), mutations.getToState(i),
                                (i - start + 1.0) / (end - start + 1.0));
                    }
                    start = end;
                }
                for (int i = mutations.getSubstitutionCount(); i < mutations.size(); i++) {
                    builder.addAmbiguity(mutations.getPattern(i), mutations.getFromState(i));
                }
                store.setBranchMutations(n, builder.build());
            }
        }
        store.commit();
    }
}
