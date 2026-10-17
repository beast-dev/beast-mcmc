/*
 * BranchHistoryModel.java
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

import dr.evomodel.substmodel.SubstitutionModel;

/**
 * How the mutational history of each site along a branch is represented and scored. The MAT
 * likelihood and the mutation history resampler delegate everything that depends on this to the model:
 *
 * <ul>
 * <li>{@link SingleHitHistoryModel}: at most one substitution per site per branch, times integrated out</li>
 * <li>{@link IntegratedHitsHistoryModel}: any number of hits, summed analytically (transition probabilities)</li>
 * <li>{@link ExplicitHitsHistoryModel}: explicit hits (paths of states) per site per branch, times integrated out</li>
 * <li>{@link TimedHitsHistoryModel}: explicit hits with their positions along the branch (as in Varilly et al., Nature 2026)</li>
 * </ul>
 *
 * Branch lengths are in expected substitutions (time x branch rate) and nu is a site's relative rate.
 * A site that becomes ambiguous on a branch contributes nothing to that branch.
 */
public interface BranchHistoryModel {

    String[] NAMES = {SingleHitHistoryModel.NAME, IntegratedHitsHistoryModel.NAME,
            ExplicitHitsHistoryModel.NAME, TimedHitsHistoryModel.NAME};

    /**
     * @return the branch history model with the given name (one of {@link #NAMES})
     */
    static BranchHistoryModel create(String name, int stateCount) {
        if (SingleHitHistoryModel.NAME.equalsIgnoreCase(name)) {
            return new SingleHitHistoryModel(stateCount);
        } else if (IntegratedHitsHistoryModel.NAME.equalsIgnoreCase(name)) {
            return new IntegratedHitsHistoryModel(stateCount);
        } else if (ExplicitHitsHistoryModel.NAME.equalsIgnoreCase(name)) {
            return new ExplicitHitsHistoryModel(stateCount);
        } else if (TimedHitsHistoryModel.NAME.equalsIgnoreCase(name)) {
            return new TimedHitsHistoryModel(stateCount);
        }
        throw new IllegalArgumentException("Unknown branch history model '" + name + "' (should be one of " +
                String.join(", ", NAMES) + ")");
    }

    String getName();

    /**
     * Called whenever the rate matrix may have changed.
     *
     * @param infinitesimalMatrix the (normalised) rate matrix Q, row-major
     * @param exitRates           q_a = sum over b != a of Q_ab
     */
    void setRates(SubstitutionModel substitutionModel, double[] infinitesimalMatrix, double[] exitRates);

    /**
     * Checks that the model can be used with the likelihood's settings.
     *
     * @throws IllegalArgumentException if not
     */
    void checkCompatibility(MATDataLikelihoodDelegate delegate);

    /**
     * @return true if within-branch histories (paths of hits) are part of the sampled state, so that the
     * resampler must draw them along with the ancestral states
     */
    boolean hasExplicitHistories();

    /**
     * @return true if branch terms need the state composition of each node (the summed weights of the
     * sites in each state)
     */
    boolean needsComposition();

    /**
     * @return true if branch terms need the exit rate of each node (lambda: the sum of w_p nu_p q(state of
     * p) over its sites)
     */
    default boolean needsExitRate() {
        return true;
    }

    /**
     * The log density of going from state a to state b along a branch, summed over any part of the
     * history that the model does not sample (used to draw ancestral states).
     *
     * @param out filled with the densities, row-major [a * stateCount + b]
     */
    void getLogTransitionDensities(double length, double nu, double[] out);

    /**
     * @param exitRateSum   lambda: the sum of w_p nu_p q(state of p) over the sites present at the node
     * @param composition   the summed weights of the sites present at the node in each state (only if
     *                      {@link #needsComposition()})
     * @return the log density of the histories of all sites present at the node on the branch above it
     */
    double getBranchLogLikelihood(BranchMutations mutations, double length, double exitRateSum, double[] composition,
                                  double[] patternWeights, double[] patternRates);

    /**
     * @return the log density of the history of a site with no events on a branch, in state a throughout
     */
    double getUnchangedSiteLogLikelihood(int a, double length, double nu);

    /**
     * @return the log density of the history of a site given by its substitutions [start, end) on a branch
     */
    double getSiteLogLikelihood(BranchMutations mutations, int start, int end, double length, double nu);

    /**
     * @return the log probability of a site's history on a branch given its end states (start == end
     * for no substitutions); 0 for models without explicit histories
     */
    double getHistoryLogProposal(BranchMutations mutations, int start, int end, int fromState, int toState,
                                 double length, double nu);

    /**
     * Draws a history for a site on a branch from fromState to toState, adding its substitutions to the
     * builder (models without explicit histories add the single substitution if the states differ).
     *
     * @return the log probability of the drawn history given its end states
     */
    double sampleHistory(int pattern, int fromState, int toState, double length, double nu, BranchMutations.Builder builder);

    /**
     * Gives every branch's events the form the model expects (e.g. positions for timed hits), after they
     * have been built by parsimony.
     */
    void initialise(DefaultMutationStore store);
}
