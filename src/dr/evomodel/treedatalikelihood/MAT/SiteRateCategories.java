/*
 * SiteRateCategories.java
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

import dr.evolution.alignment.PatternList;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.math.MathUtils;

import java.util.Arrays;
import java.util.logging.Logger;

/**
 * The discrete site rate categories of a MAT likelihood whose site rate model has more than one category
 * (see optimization_notes/mat_site_rate_categories/DESIGN.md).
 * <p>
 * Each variable pattern's category is part of the sampled state (kept by the store, journaled, and drawn
 * by {@link #resample}, a Gibbs update), so that the variable sites' terms still decompose by branch, each
 * at its category's rate; summed over the categories the likelihood is that of the site rate model. The
 * constant patterns are summed over the categories analytically: each has no events, so its likelihood
 * in a category is exp of the sum over the branches of the unchanged-site term, which this class keeps
 * as running totals over the branches.
 */
final class SiteRateCategories {

    private final MATDataLikelihoodDelegate delegate;
    private final SiteRateModel siteRateModel;
    private final Tree tree;
    private final DefaultMutationStore store;
    private final BranchHistoryModel historyModel;
    private final int categoryCount;
    private final int stateCount;
    private final int patternCount;
    private final double[] patternWeights;

    private final double[] rates;                  // the rate of each category (see updateRates)
    private final double[] logProportions;         // the log prior probability of each category
    private final double[] constantStateWeights;   // the summed weights of the constant patterns in each state

    // the unchanged-site terms of the branch above each node, [node][category * stateCount + state], and
    // their running totals over the branches (each node's terms are included at most once)
    private final double[][] unchangedTerms;
    private final double[] unchangedTermLengths;
    private final long[] unchangedTermModelVersions;
    private final boolean[] unchangedTermIncluded;
    private final double[] unchangedTotals;
    private int unchangedExcludedNode = -1;        // the root when its term was taken out of the totals

    private double priorTerm;
    private long priorModelVersion = -1;
    private boolean hasOperator = false;
    private boolean warningGiven = false;

    /**
     * Gives every variable pattern of the delegate's store the category whose rate is closest to 1.
     */
    SiteRateCategories(MATDataLikelihoodDelegate delegate, SiteRateModel siteRateModel, PatternList patternList) {
        this.delegate = delegate;
        this.siteRateModel = siteRateModel;
        this.tree = delegate.getTree();
        this.store = delegate.getMutableMutationStore();
        this.historyModel = delegate.getHistoryModel();
        this.categoryCount = siteRateModel.getCategoryCount();
        this.stateCount = patternList.getDataType().getStateCount();
        this.patternCount = patternList.getPatternCount();
        this.patternWeights = delegate.getPatternWeights();
        this.rates = new double[categoryCount];
        this.logProportions = new double[categoryCount];
        this.constantStateWeights = new double[stateCount];

        int initial = -1;
        for (int c = 0; c < categoryCount; c++) {
            double rate = siteRateModel.getRateForCategory(c);
            if (rate > 0.0 && (initial < 0 ||
                    Math.abs(Math.log(rate)) < Math.abs(Math.log(siteRateModel.getRateForCategory(initial))))) {
                initial = c;
            }
        }
        byte[] categories = new byte[patternCount];
        for (int p = 0; p < patternCount; p++) {
            boolean constant = delegate.isConstantPattern(p);
            categories[p] = (byte) (constant ? -1 : initial);
            if (constant) {
                constantStateWeights[patternList.getPatternState(0, p)] += patternWeights[p];
            }
        }
        store.setCategories(categoryCount, categories);

        final int nodeCount = tree.getNodeCount();
        this.unchangedTerms = new double[nodeCount][categoryCount * stateCount];
        this.unchangedTermLengths = new double[nodeCount];
        this.unchangedTermModelVersions = new long[nodeCount];
        Arrays.fill(unchangedTermModelVersions, -1);
        this.unchangedTermIncluded = new boolean[nodeCount];
        this.unchangedTotals = new double[categoryCount * stateCount];
        Logger.getLogger("dr.evomodel").info("  MAT site rate categories: " + categoryCount +
                " (sampled for each variable site; constant sites summed over them)");
    }

    int getCategoryCount() {
        return categoryCount;
    }

    /**
     * @return the rate of each category (after {@link #updateRates}); callers must not modify it
     */
    double[] getRates() {
        return rates;
    }

    /**
     * Reads the categories' rates and proportions from the site rate model.
     *
     * @return true if any category's rate changed
     */
    boolean updateRates() {
        boolean changed = false;
        for (int c = 0; c < categoryCount; c++) {
            double rate = siteRateModel.getRateForCategory(c);
            if (rate != rates[c]) {
                rates[c] = rate;
                changed = true;
            }
            logProportions[c] = Math.log(siteRateModel.getProportionForCategory(c));
        }
        return changed;
    }

    /**
     * @return the rate of the pattern's category, or 0 for a constant pattern (which is summed over the
     * categories separately, see {@link #getConstantTerm})
     */
    double getPatternRate(int pattern) {
        return delegate.isConstantPattern(pattern) ? 0.0 : rates[store.getCategory(pattern)];
    }

    /**
     * Called by a {@link SiteRateCategoryGibbsOperator} that updates the categories.
     */
    void setOperator() {
        hasOperator = true;
    }

    void warnIfNoOperator() {
        if (!hasOperator && !warningGiven) {
            // (at the info level, as BEAST prints warnings only when verbose)
            Logger.getLogger("dr.evomodel").info("WARNING: the MAT likelihood samples a site rate category for " +
                    "each variable site, but no siteRateCategoryGibbsOperator updates them: add one to the operators " +
                    "(with a weight of about 1% of the subtree moves' weight).");
            warningGiven = true;
        }
    }

    /**
     * Empties the running totals, to be rebuilt from every node's terms (which also stops rounding errors
     * accumulating).
     */
    void resetTotals() {
        Arrays.fill(unchangedTotals, 0.0);
        Arrays.fill(unchangedTermIncluded, false);
        unchangedExcludedNode = -1;
    }

    /**
     * Brings the node's unchanged-site terms (for the branch above it) up to date, and includes them in the
     * totals.
     *
     * @param ratesVersion changes whenever the rate matrix or the categories' rates change (the terms do not
     *                     depend on the states or the categories)
     */
    void updateUnchangedTerm(int n, double length, long ratesVersion) {
        final double[] terms = unchangedTerms[n];
        final boolean stale = unchangedTermModelVersions[n] != ratesVersion || unchangedTermLengths[n] != length;
        if (unchangedTermIncluded[n] && !stale) {
            return;
        }
        if (unchangedTermIncluded[n]) {
            for (int i = 0; i < terms.length; i++) {
                unchangedTotals[i] -= terms[i];
            }
        }
        if (stale) {
            for (int c = 0; c < categoryCount; c++) {
                historyModel.getUnchangedSiteLogLikelihoods(length, rates[c], terms, c * stateCount);
            }
            unchangedTermLengths[n] = length;
            unchangedTermModelVersions[n] = ratesVersion;
        }
        for (int i = 0; i < terms.length; i++) {
            unchangedTotals[i] += terms[i];
        }
        unchangedTermIncluded[n] = true;
    }

    /**
     * @return the node's unchanged-site terms, [category * stateCount + state] (after {@link #updateUnchangedTerm})
     */
    double[] getUnchangedTerms(int n) {
        return unchangedTerms[n];
    }

    /**
     * @return the log likelihood of the constant patterns, summed over the site rate categories: each has no
     * events, so in category c and state a it is the sum of the unchanged-site terms over the branches
     */
    double getConstantTerm(int rootNumber, long ratesVersion) {
        // the totals cover every branch: the root's (former) branch is taken out, and a former root's put back
        if (unchangedExcludedNode >= 0 && unchangedExcludedNode != rootNumber) {
            int n = unchangedExcludedNode;
            unchangedExcludedNode = -1;
            updateUnchangedTerm(n, delegate.getBranchLength(tree.getNode(n)), ratesVersion);
        }
        if (unchangedTermIncluded[rootNumber]) {
            final double[] terms = unchangedTerms[rootNumber];
            for (int i = 0; i < terms.length; i++) {
                unchangedTotals[i] -= terms[i];
            }
            unchangedTermIncluded[rootNumber] = false;
        }
        unchangedExcludedNode = rootNumber;

        double term = 0.0;
        for (int a = 0; a < stateCount; a++) {
            if (constantStateWeights[a] > 0.0) {
                term += constantStateWeights[a] * logSumCategories(unchangedTotals, a);
            }
        }
        return term;
    }

    /**
     * @return log of the sum over categories of their proportion times exp(values[c * stateCount + a])
     */
    private double logSumCategories(double[] values, int a) {
        double max = Double.NEGATIVE_INFINITY;
        for (int c = 0; c < categoryCount; c++) {
            max = Math.max(max, logProportions[c] + values[c * stateCount + a]);
        }
        if (max == Double.NEGATIVE_INFINITY) {
            return max;
        }
        double sum = 0.0;
        for (int c = 0; c < categoryCount; c++) {
            sum += Math.exp(logProportions[c] + values[c * stateCount + a] - max);
        }
        return max + Math.log(sum);
    }

    /**
     * @return the log prior probability of the variable patterns' categories
     */
    double getPriorTerm(long modelVersion) {
        if (priorModelVersion != modelVersion) {
            int[] counts = new int[categoryCount];
            for (int p = 0; p < patternCount; p++) {
                if (!delegate.isConstantPattern(p)) {
                    counts[store.getCategory(p)]++;
                }
            }
            double term = 0.0;
            for (int c = 0; c < categoryCount; c++) {
                if (counts[c] > 0) {
                    term += counts[c] * logProportions[c];
                }
            }
            priorTerm = term;
            priorModelVersion = modelVersion;
        }
        return priorTerm;
    }

    /**
     * Completes the per-pattern log likelihoods (computed in full by the delegate, each variable pattern at
     * its category's rate): adds the variable patterns' category prior, and sums the constant patterns over
     * the categories.
     */
    void completeSiteLogLikelihoods(double[] siteLogLikelihoods, double[] frequencies, byte[] rootStates,
                                    int rootNumber) {
        double[] totals = new double[categoryCount * stateCount];
        for (int n = 0; n < tree.getNodeCount(); n++) {
            if (n != rootNumber) {
                double length = delegate.getBranchLength(tree.getNode(n));
                double[] terms = new double[categoryCount * stateCount];
                for (int c = 0; c < categoryCount; c++) {
                    historyModel.getUnchangedSiteLogLikelihoods(length, rates[c], terms, c * stateCount);
                }
                for (int i = 0; i < terms.length; i++) {
                    totals[i] += terms[i];
                }
            }
        }
        for (int p = 0; p < patternCount; p++) {
            if (delegate.isConstantPattern(p)) {
                siteLogLikelihoods[p] = Math.log(frequencies[rootStates[p]]) + logSumCategories(totals, rootStates[p]);
            } else {
                siteLogLikelihoods[p] += logProportions[store.getCategory(p)];
            }
        }
    }

    /**
     * Draws the site rate category of every variable pattern from its full conditional given the history
     * and the tree (a Gibbs update): in proportion to the category's prior probability times the site's
     * likelihood at its rate. The MAT and the rates must be up to date.
     * <p>
     * A site's log likelihood at rate r is that of its root state staying unchanged on every branch, plus a
     * correction at each of its events, since an event on the branch above node n changes the state of
     * the whole subtree below n (until the next event). So with D_s(n; r), the sum of the unchanged-site
     * terms in state s over the branches below n, the likelihoods of every site at every rate take
     * O(nodes x categories x states + events x categories).
     */
    void resample() {
        final int nodeCount = tree.getNodeCount();
        final int width = categoryCount * stateCount;
        final double[][] unchanged = new double[nodeCount][width]; // for the branch above each node
        final double[][] below = new double[nodeCount][width];     // D_s(n; r_c), summed below each node
        final double[] lengths = new double[nodeCount];

        final int[] order = MATInitializer.getPostOrder(tree);
        final int rootNumber = tree.getRoot().getNumber();
        for (int n : order) {
            NodeRef node = tree.getNode(n);
            for (int i = 0; i < tree.getChildCount(node); i++) {
                int child = tree.getChild(node, i).getNumber();
                for (int j = 0; j < width; j++) {
                    below[n][j] += unchanged[child][j] + below[child][j];
                }
            }
            if (n != rootNumber) {
                lengths[n] = delegate.getBranchLength(node);
                for (int c = 0; c < categoryCount; c++) {
                    historyModel.getUnchangedSiteLogLikelihoods(lengths[n], rates[c], unchanged[n], c * stateCount);
                }
            }
        }

        // each site's log likelihood (less its root term) in each category
        final double[] logL = new double[patternCount * categoryCount];
        final byte[] rootStates = store.getStates(rootNumber);
        for (int p = 0; p < patternCount; p++) {
            if (!delegate.isConstantPattern(p) && rootStates[p] >= 0) {
                for (int c = 0; c < categoryCount; c++) {
                    logL[p * categoryCount + c] = below[rootNumber][c * stateCount + rootStates[p]];
                }
            }
        }
        for (int n = 0; n < nodeCount; n++) {
            if (n == rootNumber) {
                continue;
            }
            final BranchMutations mutations = store.getBranchMutations(n);
            final double[] u = unchanged[n];
            final double[] d = below[n];
            int i = 0;
            while (i < mutations.getSubstitutionCount()) {
                // the site's state changes from a to b on this branch (and below)
                final int end = mutations.getSiteEnd(i);
                final int p = mutations.getPattern(i);
                final int a = mutations.getFromState(i);
                final int b = mutations.getToState(end - 1);
                for (int c = 0; c < categoryCount; c++) {
                    double events = historyModel.getSiteLogLikelihood(mutations, i, end, lengths[n], rates[c]);
                    logL[p * categoryCount + c] += events - u[c * stateCount + a]
                            - d[c * stateCount + a] + d[c * stateCount + b];
                }
                i = end;
            }
            for (; i < mutations.size(); i++) {
                // no data below: the branch and the subtree contribute nothing
                final int p = mutations.getPattern(i);
                final int a = mutations.getFromState(i);
                for (int c = 0; c < categoryCount; c++) {
                    logL[p * categoryCount + c] -= u[c * stateCount + a] + d[c * stateCount + a];
                }
            }
        }

        final double[] weights = new double[categoryCount];
        final byte[] categories = store.getCategories().clone();
        for (int p = 0; p < patternCount; p++) {
            if (delegate.isConstantPattern(p)) {
                continue;
            }
            double max = Double.NEGATIVE_INFINITY;
            for (int c = 0; c < categoryCount; c++) {
                weights[c] = logProportions[c] + patternWeights[p] * logL[p * categoryCount + c];
                max = Math.max(max, weights[c]);
            }
            double sum = 0.0;
            for (int c = 0; c < categoryCount; c++) {
                weights[c] = Math.exp(weights[c] - max);
                sum += weights[c];
            }
            double u = MathUtils.nextDouble() * sum;
            int chosen = categoryCount - 1;
            for (int c = 0; c < categoryCount; c++) {
                u -= weights[c];
                if (u < 0.0) {
                    chosen = c;
                    break;
                }
            }
            categories[p] = (byte) chosen;
        }

        // the nodes' category compositions for the new categories, all at once: the root's from its states,
        // each other node's from its parent's and the events on its branch (in pre-order)
        final double[][] compositions = new double[nodeCount][];
        compositions[rootNumber] = new double[width];
        for (int p = 0; p < patternCount; p++) {
            if (!delegate.isConstantPattern(p) && rootStates[p] >= 0) {
                compositions[rootNumber][categories[p] * stateCount + rootStates[p]] += patternWeights[p];
            }
        }
        for (int j = nodeCount - 1; j >= 0; j--) {
            final int n = order[j];
            if (n == rootNumber) {
                continue;
            }
            final double[] composition = compositions[store.getParent(n)].clone();
            final BranchMutations mutations = store.getBranchMutations(n);
            int i = 0;
            while (i < mutations.getSubstitutionCount()) {
                final int end = mutations.getSiteEnd(i);
                final int p = mutations.getPattern(i);
                final int offset = categories[p] * stateCount;
                composition[offset + mutations.getFromState(i)] -= patternWeights[p];
                composition[offset + mutations.getToState(end - 1)] += patternWeights[p];
                i = end;
            }
            for (; i < mutations.size(); i++) {
                final int p = mutations.getPattern(i);
                composition[categories[p] * stateCount + mutations.getFromState(i)] -= patternWeights[p];
            }
            compositions[n] = composition;
        }
        store.replaceCategories(categories, compositions);
    }
}
