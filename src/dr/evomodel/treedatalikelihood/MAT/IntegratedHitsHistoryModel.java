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

import dr.evomodel.substmodel.BaseSubstitutionModel;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.SubstitutionModel;

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
    public void setRates(SubstitutionModel substitutionModel, double[] infinitesimalMatrix, double[] exitRates) {
        super.setRates(substitutionModel, infinitesimalMatrix, exitRates);
        setEigenDecomposition(substitutionModel);
    }

    @Override
    public void checkCompatibility(MATDataLikelihoodDelegate delegate) {
        if (delegate.hasSiteRates()) {
            throw new IllegalArgumentException("The " + NAME + " branch history model does not support per-site rates");
        }
    }

    @Override
    public boolean needsComposition() {
        return true;
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
    public double getCategoryBranchLogLikelihood(BranchMutations mutations, double length, double exitRateSum,
                                                 double[] categoryComposition, double[] unchangedTerms,
                                                 double[] patternWeights, double[] patternRates) {
        // the unchanged sites: the diagonals log P_aa(r_c t)
        double term = 0.0;
        for (int i = 0; i < categoryComposition.length; i++) {
            if (categoryComposition[i] > 0.0) {
                term += categoryComposition[i] * unchangedTerms[i];
            }
        }
        final boolean direct = hasDirectTransitionProbabilities();
        for (int k = 0; k < mutations.getSubstitutionCount(); k++) {
            int p = mutations.getPattern(k);
            int a = mutations.getFromState(k);
            int b = mutations.getToState(k);
            double distance = patternRates[p] * length;
            // the composition counted the site as staying in state b
            if (direct) {
                term += patternWeights[p] * (getLogTransitionProbability(distance, a, b) -
                        getLogTransitionProbability(distance, b, b));
            } else {
                double[] logP = getLogTransitionProbabilities(distance);
                term += patternWeights[p] * (logP[a * stateCount + b] - logP[b * stateCount + b]);
            }
        }
        return term;
    }

    @Override
    public void getUnchangedSiteLogLikelihoods(double length, double nu, double[] out, int offset) {
        if (hasDirectTransitionProbabilities()) {
            // only the diagonal (with the cache's values)
            getLogTransitionDiagonal(nu * length, out, offset);
            return;
        }
        // one transition matrix for all the states
        double[] logP = getLogTransitionProbabilities(nu * length);
        for (int a = 0; a < stateCount; a++) {
            out[offset + a] = logP[a * stateCount + a];
        }
    }

    @Override
    public double getUnchangedSiteLogLikelihood(int a, double length, double nu) {
        return getLogTransitionProbabilities(nu * length)[a * stateCount + a];
    }

    @Override
    public double getSiteLogLikelihood(BranchMutations mutations, int start, int end, double length, double nu) {
        return getLogTransitionProbabilities(nu * length)[mutations.getFromState(start) * stateCount + mutations.getToState(end - 1)];
    }

    // ---------------------------------------------------------------------------------------------
    // Transition probabilities from the eigen decomposition (without the cache), for the site rate
    // categories' terms
    // ---------------------------------------------------------------------------------------------

    // the rate matrix's eigen decomposition, when the substitution model computes its transition
    // probabilities from it as BaseSubstitutionModel does (otherwise null)
    private double[] eigenValues;
    private double[] eigenVectors;
    private double[] inverseEigenVectors;
    private double[] eigenExp;

    private void setEigenDecomposition(SubstitutionModel substitutionModel) {
        eigenValues = null;
        if (!(substitutionModel instanceof BaseSubstitutionModel)) {
            return;
        }
        try {
            if (substitutionModel.getClass().getMethod("getTransitionProbabilities", double.class, double[].class)
                    .getDeclaringClass() != BaseSubstitutionModel.class) {
                return;
            }
        } catch (NoSuchMethodException e) {
            return;
        }
        EigenDecomposition eigen = substitutionModel.getEigenDecomposition();
        if (eigen == null || eigen.getEigenValues().length != stateCount ||
                eigen.getEigenVectors().length != stateCount * stateCount) {
            return;
        }
        eigenValues = eigen.getEigenValues().clone();
        eigenVectors = eigen.getEigenVectors().clone();
        inverseEigenVectors = eigen.getInverseEigenVectors().clone();
        eigenExp = new double[stateCount];
    }

    /**
     * @return true if transition probabilities can be computed directly ({@link #getLogTransitionDiagonal},
     * {@link #getLogTransitionProbability})
     */
    private boolean hasDirectTransitionProbabilities() {
        return eigenValues != null;
    }

    /**
     * Sets out[offset + a] to log P_aa(distance) for every state a: the same values as the diagonal of
     * {@link #getLogTransitionProbabilities} (computed with the same arithmetic as BaseSubstitutionModel),
     * for S exponentials and logarithms rather than a whole matrix.
     */
    private void getLogTransitionDiagonal(double distance, double[] out, int offset) {
        for (int i = 0; i < stateCount; i++) {
            eigenExp[i] = Math.exp(distance * eigenValues[i]);
        }
        for (int a = 0; a < stateCount; a++) {
            double p = 0.0;
            for (int k = 0; k < stateCount; k++) {
                p += eigenVectors[a * stateCount + k] * (inverseEigenVectors[k * stateCount + a] * eigenExp[k]);
            }
            out[offset + a] = Math.log(Math.abs(p));
        }
    }

    /**
     * @return log P_ab(distance), as {@link #getLogTransitionProbabilities} gives it
     */
    private double getLogTransitionProbability(double distance, int a, int b) {
        double p = 0.0;
        for (int k = 0; k < stateCount; k++) {
            p += eigenVectors[a * stateCount + k] * (inverseEigenVectors[k * stateCount + b] * Math.exp(distance * eigenValues[k]));
        }
        return Math.log(Math.abs(p));
    }
}
