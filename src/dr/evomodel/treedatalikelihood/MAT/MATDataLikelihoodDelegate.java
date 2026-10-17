/*
 * MATDataLikelihoodDelegate.java
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
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.RateRescalingScheme;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.TreeTraversal;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.List;
import java.util.logging.Logger;

/**
 * Mutation-annotated tree (MAT) likelihood delegate, following the model used by Delphy
 * (Varilly et al.).
 *
 * Rather than integrating over all ancestral states, the likelihood is that of an explicit
 * mutational history: a root sequence plus a set of mutations on each branch, reconstructed by
 * parsimony ({@link MATInitializer}). The ancestral states can then be sampled jointly with the
 * topology by SPR operators carrying a {@link MutationHistoryResampler}. If the topology is changed by
 * any other operator, the MAT is brought up to date by (incremental) parsimony, which makes the
 * likelihood a deterministic function of the tree. The two should not be mixed in one analysis.
 *
 * How each site's history along a branch is represented and scored is delegated to a
 * {@link BranchHistoryModel}: by default {@link SingleHitHistoryModel} (at most one substitution per
 * site per branch, its time integrated out); alternatives allow multiple hits, summed analytically or
 * as explicit (optionally timed) paths. The root sequence has prior sum_l log pi_a.
 *
 * The likelihood is a root term plus a term for each branch, each cached against what it depends on,
 * so an evaluation costs O(changed branches) plus O(nodes) to check the caches.
 */
public class MATDataLikelihoodDelegate extends AbstractModel implements DataLikelihoodDelegate {

    public static final String NAME = "MATDataLikelihoodDelegate";

    private final PatternList patternList;
    private final Tree tree;
    private final BranchModel branchModel;
    private final SubstitutionModel substitutionModel;
    private final FrequencyModel rootFrequencyModel;
    private final SiteRateModel siteRateModel;
    private final BranchRateModel branchRateModel;
    private final Parameter siteRates;

    private final int patternCount;
    private final int stateCount;
    private final double[] patternWeights;
    private final boolean[] constantPatterns;

    private final BranchHistoryModel historyModel;
    private final MATInitializer initializer;
    private final DefaultMutationStore mutationStore;
    private int rebuildCount = 0;
    private boolean historyResampled = false;
    private boolean mixingWarningGiven = false;

    private final double[] siteLogLikelihoods;
    private boolean siteLogLikelihoodsKnown = false;
    private double logLikelihood;
    private double storedLogLikelihood;

    // Incremental likelihood: logL = root term + sum over branches of a branch term, each cached
    // against what it depends on. modelVersion changes whenever the rate matrix, frequencies or
    // site rates may have changed (branch lengths are compared directly).
    private long modelVersion = 0;
    private long cachedModelVersion = -1;
    private boolean modelChangedSinceStore = false;
    private final double[] patternRates;        // nu_p
    private final double[] patternWeightRates;  // w_p nu_p
    private final double[] exitRateSums;        // lambda_n = sum_p w_p nu_p q(state of p at n)
    private final long[] exitRateStateVersions;
    private final int[] exitRatePath;           // scratch: nodes whose exit rates are being brought up to date
    private final double[] branchTerms;
    private final double[] branchTermLengths;
    private final long[] branchTermStateVersions;
    private double commonSiteRate; // the rate of every site, without per-site rates
    private final BranchMutations[] branchTermMutations;
    private final long[] branchTermModelVersions;
    private final long[] exitRateModelVersions;
    private double rootTerm;
    private int rootTermNode = -1;
    private long rootTermStateVersion = -1;
    private long rootTermModelVersion = -1;
    private long calculationCount = 0;

    // The branch terms are only recomputed for nodes whose states, branch events or branch length may
    // have changed: the store logs the first two, events from the tree the third.
    private final boolean[] nodeDirty;
    private final int[] dirtyNodes;
    private int dirtyNodeCount = 0;
    private boolean allNodesDirty = true;
    private long termsModelVersion = -1;
    // nodes whose terms were recomputed since storeState (the tree restores their lengths silently)
    private final boolean[] recomputedSinceStore;
    private final int[] recomputedNodes;
    private int recomputedCount = 0;
    // whether the topology may have changed since the MAT was last checked against it
    private boolean treeChangedSinceCheck = true;

    private final double[] infinitesimalMatrix;
    private final double[] exitRates;

    public MATDataLikelihoodDelegate(Tree tree,
                                     PatternList patternList,
                                     BranchModel branchModel,
                                     SiteRateModel siteRateModel,
                                     BranchRateModel branchRateModel,
                                     Parameter siteRates) {
        this(tree, patternList, branchModel, siteRateModel, branchRateModel, siteRates, null);
    }

    /**
     * @param historyModel how histories along branches are represented (null for single hits)
     */
    public MATDataLikelihoodDelegate(Tree tree,
                                     PatternList patternList,
                                     BranchModel branchModel,
                                     SiteRateModel siteRateModel,
                                     BranchRateModel branchRateModel,
                                     Parameter siteRates,
                                     BranchHistoryModel historyModel) {
        super(NAME);
        this.tree = tree;
        this.patternList = patternList;
        this.branchModel = branchModel;
        this.siteRateModel = siteRateModel;
        this.branchRateModel = branchRateModel == null ? new DefaultBranchRateModel() : branchRateModel;
        this.siteRates = siteRates;

        if (branchModel.getSubstitutionModels().size() != 1) {
            throw new IllegalArgumentException("MAT likelihood currently requires a single substitution model across the tree");
        }
        this.substitutionModel = branchModel.getSubstitutionModels().get(0);
        this.rootFrequencyModel = branchModel.getRootFrequencyModel();

        if (siteRateModel != null && siteRateModel.getCategoryCount() != 1) {
            throw new IllegalArgumentException("MAT likelihood does not support discrete site rate categories; " +
                    "use a per-site rate parameter instead");
        }

        this.patternCount = patternList.getPatternCount();
        this.stateCount = patternList.getDataType().getStateCount();
        if (substitutionModel.getDataType().getStateCount() != stateCount) {
            throw new IllegalArgumentException("Substitution model and pattern list have different state counts");
        }
        if (siteRates != null && siteRates.getDimension() != patternCount) {
            throw new IllegalArgumentException("Site rate parameter " + siteRates.getId() + " has dimension " +
                    siteRates.getDimension() + " but there are " + patternCount + " patterns");
        }

        this.patternWeights = patternList.getPatternWeights();

        // Each site has its own mutational history, so the patterns must not be compressed, except that
        // constant sites (which are assumed to have no mutations) may be pooled by state, as produced by
        // the constantOnly compression of SitePatterns.
        this.constantPatterns = new boolean[patternCount];
        int constantSites = 0;
        int variableSites = 0;
        for (int p = 0; p < patternCount; p++) {
            int first = patternList.getPatternState(0, p);
            boolean constant = first < stateCount;
            for (int i = 1; constant && i < patternList.getPatternLength(); i++) {
                constant = patternList.getPatternState(i, p) == first;
            }
            constantPatterns[p] = constant;
            if (constant) {
                constantSites += (int) patternWeights[p];
            } else {
                if (patternWeights[p] != 1.0) {
                    throw new IllegalArgumentException("MAT likelihood needs each variable site as its own pattern but pattern " +
                            (p + 1) + " has weight " + patternWeights[p] + ": use uncompressed patterns or constant-only " +
                            "compression (unique=\"false\" or constantOnly=\"true\" on the patterns element)");
                }
                variableSites++;
            }
        }
        Logger.getLogger("dr.evomodel").info("  MAT likelihood: " + variableSites + " variable sites, " + constantSites +
                " constant sites (assumed to have no mutations)");
        this.siteLogLikelihoods = new double[patternCount];
        int nodeCount = tree.getNodeCount();
        this.patternRates = new double[patternCount];
        this.patternWeightRates = new double[patternCount];
        this.exitRateSums = new double[nodeCount];
        this.exitRateStateVersions = new long[nodeCount];
        this.exitRateModelVersions = new long[nodeCount];
        this.exitRatePath = new int[nodeCount];
        this.branchTerms = new double[nodeCount];
        this.branchTermLengths = new double[nodeCount];
        this.branchTermStateVersions = new long[nodeCount];
        this.branchTermMutations = new BranchMutations[nodeCount];
        this.branchTermModelVersions = new long[nodeCount];
        java.util.Arrays.fill(exitRateModelVersions, -1);
        java.util.Arrays.fill(branchTermModelVersions, -1);
        this.nodeDirty = new boolean[nodeCount];
        this.dirtyNodes = new int[nodeCount];
        this.recomputedSinceStore = new boolean[nodeCount];
        this.recomputedNodes = new int[nodeCount];
        this.infinitesimalMatrix = new double[stateCount * stateCount];
        this.exitRates = new double[stateCount];

        this.historyModel = historyModel == null ? new SingleHitHistoryModel(stateCount) : historyModel;
        this.historyModel.checkCompatibility(this);
        Logger.getLogger("dr.evomodel").info("  MAT branch history model: " + this.historyModel.getName());

        initializer = new MATInitializer(patternList, tree);
        mutationStore = initializer.buildMutationStore();
        this.historyModel.initialise(mutationStore);

        addModel(branchModel);
        if (siteRateModel != null) {
            addModel(siteRateModel);
        }
        addModel(this.branchRateModel);
        if (siteRates != null) {
            addVariable(siteRates);
        }
        if (tree instanceof Model) {
            // a listener only (not a sub-model): to learn which branch lengths, and when the topology, change
            ((Model) tree).addModelListener(this);
        }
    }

    public MutationStore getMutationStore() {
        return mutationStore;
    }

    public BranchHistoryModel getHistoryModel() {
        return historyModel;
    }

    MATInitializer getInitializer() {
        return initializer;
    }

    DefaultMutationStore getMutableMutationStore() {
        return mutationStore;
    }

    /**
     * Called by a {@link MutationHistoryResampler} when it has sampled new ancestral states.
     */
    void setHistoryResampled() {
        historyResampled = true;
    }

    /**
     * @return the number of times the MAT has been updated by parsimony after a topology change
     */
    public int getRebuildCount() {
        return rebuildCount;
    }

    public Tree getTree() {
        return tree;
    }

    public PatternList getPatternList() {
        return patternList;
    }

    public BranchModel getBranchModel() {
        return branchModel;
    }

    public SiteRateModel getSiteRateModel() {
        return siteRateModel;
    }

    public BranchRateModel getBranchRateModel() {
        return branchRateModel;
    }

    public Parameter getSiteRates() {
        return siteRates;
    }

    /**
     * Makes the MAT consistent with the current topology, by parsimony if necessary.
     */
    void updateMAT() {
        if (!treeChangedSinceCheck) {
            return;
        }
        treeChangedSinceCheck = false;
        if (!mutationStore.hasParents(MATInitializer.getParents(tree))) {
            initializer.update(mutationStore);
            historyModel.initialise(mutationStore);
            rebuildCount++;
            if (historyResampled && !mixingWarningGiven) {
                Logger.getLogger("dr.evomodel").warning("The topology of " + tree.getId() +
                        " was changed by an operator without a mutation history resampler: the MAT has been " +
                        "updated by parsimony, which is not a valid MCMC move when histories are sampled.");
                mixingWarningGiven = true;
            }
        }
    }

    /**
     * Called by a subtree move's hook once it has made the MAT match the topology after the move, so that
     * the topology need not be checked against the MAT (it is under assertions).
     */
    void markTopologySynchronised() {
        assert mutationStore.hasParents(MATInitializer.getParents(tree)) : "the MAT does not match the topology";
        treeChangedSinceCheck = false;
    }

    /**
     * The MAT tracks changes to the tree itself (see updateMAT), so the tree likelihood need not traverse
     * the tree to collect operations.
     */
    @Override
    public boolean requiresTreeTraversal() {
        return false;
    }

    @Override
    public double calculateLikelihood(List<BranchOperation> branchOperations,
                                      List<NodeOperation> nodeOperations,
                                      int rootNodeNumber) throws LikelihoodException {
        updateMAT();
        computeLikelihood();
        calculationCount++;
        return logLikelihood;
    }

    /**
     * Brings the rate matrix, the site rates and the history model's rates up to date.
     */
    void ensureRates() {
        if (cachedModelVersion != modelVersion) {
            substitutionModel.getInfinitesimalMatrix(infinitesimalMatrix);
            for (int a = 0; a < stateCount; a++) {
                double q = 0.0;
                for (int b = 0; b < stateCount; b++) {
                    if (b != a) {
                        q += infinitesimalMatrix[a * stateCount + b];
                    }
                }
                exitRates[a] = q;
            }
            double siteRateScale = getSiteRateScale();
            commonSiteRate = siteRateScale;
            for (int p = 0; p < patternCount; p++) {
                patternRates[p] = getSiteRate(p, siteRateScale);
                patternWeightRates[p] = patternWeights[p] * patternRates[p];
            }
            historyModel.setRates(substitutionModel, infinitesimalMatrix, exitRates);
            cachedModelVersion = modelVersion;
        }
    }

    private void computeLikelihood() {
        ensureRates();
        if (termsModelVersion != modelVersion) {
            allNodesDirty = true;
            termsModelVersion = modelVersion;
        }
        for (int i = 0; i < mutationStore.getChangedNodeCount(); i++) {
            markDirty(mutationStore.getChangedNode(i));
        }
        mutationStore.clearChangeLog();

        final int rootNumber = tree.getRoot().getNumber();
        if (allNodesDirty) {
            for (int n = 0; n < tree.getNodeCount(); n++) {
                if (n != rootNumber) {
                    updateBranchTerm(n);
                }
                nodeDirty[n] = false;
            }
            allNodesDirty = false;
        } else {
            for (int i = 0; i < dirtyNodeCount; i++) {
                int n = dirtyNodes[i];
                if (n != rootNumber) {
                    updateBranchTerm(n);
                }
                nodeDirty[n] = false;
            }
        }
        dirtyNodeCount = 0;

        double total = getRootTerm(rootNumber);
        for (int n = 0; n < branchTerms.length; n++) {
            if (n != rootNumber) {
                total += branchTerms[n];
            }
        }

        logLikelihood = total;
        siteLogLikelihoodsKnown = false;
    }

    private void markDirty(int n) {
        if (!nodeDirty[n]) {
            nodeDirty[n] = true;
            dirtyNodes[dirtyNodeCount++] = n;
        }
    }

    /**
     * Recomputes the term of the branch above node n if anything it depends on has changed.
     */
    private void updateBranchTerm(int n) {
        double exitRate = historyModel.needsExitRate() ? getExitRate(n) : 0.0;
        double[] composition = historyModel.needsComposition() ? mutationStore.getComposition(n) : null;
        double length = getBranchLength(tree.getNode(n));
        BranchMutations mutations = mutationStore.getBranchMutations(n);
        long stateVersion = mutationStore.getStateVersion(n);
        if (branchTermModelVersions[n] != modelVersion || branchTermLengths[n] != length ||
                branchTermStateVersions[n] != stateVersion || branchTermMutations[n] != mutations) {
            branchTerms[n] = historyModel.getBranchLogLikelihood(mutations, length, exitRate, composition,
                    patternWeights, patternRates);
            branchTermLengths[n] = length;
            branchTermStateVersions[n] = stateVersion;
            branchTermMutations[n] = mutations;
            branchTermModelVersions[n] = modelVersion;
        }
        if (!recomputedSinceStore[n]) {
            recomputedSinceStore[n] = true;
            recomputedNodes[recomputedCount++] = n;
        }
    }

    /**
     * @return lambda_n, the total exit rate of the node's sequence: from the node's composition if all
     * sites have the same rate, otherwise from its parent's and the events on its branch (bringing the
     * rates of its ancestors up to date first, from the nearest one that is)
     */
    private double getExitRate(int n) {
        if (isExitRateKnown(n)) {
            return exitRateSums[n];
        }
        if (!hasSiteRates()) {
            double[] composition = mutationStore.getComposition(n);
            double rate = 0.0;
            for (int a = 0; a < stateCount; a++) {
                rate += composition[a] * exitRates[a];
            }
            setExitRate(n, rate * commonSiteRate);
            return exitRateSums[n];
        }
        int count = 0;
        int top = n;
        exitRatePath[count++] = top;
        while (mutationStore.getParent(top) >= 0 && !isExitRateKnown(mutationStore.getParent(top))) {
            top = mutationStore.getParent(top);
            exitRatePath[count++] = top;
        }
        for (int i = count - 1; i >= 0; i--) {
            int node = exitRatePath[i];
            int parent = mutationStore.getParent(node);
            double rate;
            if (parent < 0) {
                rate = computeExitRate(node);
            } else {
                rate = exitRateSums[parent];
                BranchMutations mutations = mutationStore.getBranchMutations(node);
                for (int k = 0; k < mutations.size(); k++) {
                    int p = mutations.getPattern(k);
                    int a = mutations.getFromState(k);
                    if (mutations.isAmbiguity(k)) {
                        rate -= patternWeightRates[p] * exitRates[a];
                    } else {
                        rate += patternWeightRates[p] * (exitRates[mutations.getToState(k)] - exitRates[a]);
                    }
                }
            }
            setExitRate(node, rate);
        }
        return exitRateSums[n];
    }

    private boolean isExitRateKnown(int n) {
        return exitRateModelVersions[n] == modelVersion && exitRateStateVersions[n] == mutationStore.getStateVersion(n);
    }

    private void setExitRate(int n, double rate) {
        exitRateSums[n] = rate;
        exitRateStateVersions[n] = mutationStore.getStateVersion(n);
        exitRateModelVersions[n] = modelVersion;
    }

    private double computeExitRate(int n) {
        byte[] states = mutationStore.getStates(n);
        double rate = 0.0;
        for (int p = 0; p < patternCount; p++) {
            if (states[p] >= 0) {
                rate += patternWeightRates[p] * exitRates[states[p]];
            }
        }
        return rate;
    }

    private double getRootTerm(int rootNumber) {
        long stateVersion = mutationStore.getStateVersion(rootNumber);
        if (rootTermNode != rootNumber || rootTermStateVersion != stateVersion || rootTermModelVersion != modelVersion) {
            double[] frequencies = rootFrequencyModel.getFrequencies();
            double[] logFrequencies = new double[stateCount];
            for (int a = 0; a < stateCount; a++) {
                logFrequencies[a] = Math.log(frequencies[a]);
            }
            double[] composition = mutationStore.getComposition(rootNumber);
            double term = 0.0;
            for (int a = 0; a < stateCount; a++) {
                if (composition[a] > 0.0) {
                    term += composition[a] * logFrequencies[a];
                }
            }
            rootTerm = term;
            rootTermNode = rootNumber;
            rootTermStateVersion = stateVersion;
            rootTermModelVersion = modelVersion;
        }
        return rootTerm;
    }

    /**
     * Computes the per-pattern log likelihoods in full (only needed when they are asked for), site by
     * site and branch by branch with the history model.
     */
    private void computeSiteLogLikelihoods() {
        ensureRates();
        double[] frequencies = rootFrequencyModel.getFrequencies();
        int rootNumber = tree.getRoot().getNumber();

        byte[] rootStates = mutationStore.getStates(rootNumber);
        for (int p = 0; p < patternCount; p++) {
            int a = rootStates[p];
            siteLogLikelihoods[p] = a < 0 ? 0.0 : Math.log(frequencies[a]);
        }
        int[] eventStart = new int[patternCount];
        int[] eventEnd = new int[patternCount];
        for (int n = 0; n < tree.getNodeCount(); n++) {
            if (n == rootNumber) {
                continue;
            }
            double length = getBranchLength(tree.getNode(n));
            byte[] states = mutationStore.getStates(n);
            BranchMutations mutations = mutationStore.getBranchMutations(n);
            java.util.Arrays.fill(eventStart, -1);
            int k = 0;
            while (k < mutations.getSubstitutionCount()) {
                int end = mutations.getSiteEnd(k);
                eventStart[mutations.getPattern(k)] = k;
                eventEnd[mutations.getPattern(k)] = end;
                k = end;
            }
            for (int p = 0; p < patternCount; p++) {
                if (states[p] < 0) {
                    continue;
                }
                siteLogLikelihoods[p] += eventStart[p] < 0 ?
                        historyModel.getUnchangedSiteLogLikelihood(states[p], length, patternRates[p]) :
                        historyModel.getSiteLogLikelihood(mutations, eventStart[p], eventEnd[p], length, patternRates[p]);
            }
        }
        for (int p = 0; p < patternCount; p++) {
            siteLogLikelihoods[p] *= patternWeights[p];
        }
        siteLogLikelihoodsKnown = true;
    }

    double[] getPatternRates() {
        ensureRates();
        return patternRates;
    }

    double getSiteRateScale() {
        return siteRateModel == null ? 1.0 : siteRateModel.getRateForCategory(0);
    }

    /**
     * @return the relative rate of the pattern, i.e. nu_l times the site rate model's rate
     */
    double getSiteRate(int pattern, double scale) {
        return siteRates == null ? scale : scale * siteRates.getParameterValue(pattern);
    }

    boolean hasSiteRates() {
        return siteRates != null;
    }

    double[] getPatternWeights() {
        return patternWeights;
    }

    /**
     * @return true if the pattern is constant (every sequence has the same unambiguous state); constant
     * patterns are assumed to have no mutations
     */
    public boolean isConstantPattern(int pattern) {
        return constantPatterns[pattern];
    }

    SubstitutionModel getSubstitutionModel() {
        return substitutionModel;
    }

    FrequencyModel getRootFrequencyModel() {
        return rootFrequencyModel;
    }

    /**
     * @return the length of the branch above the node in expected substitutions (0 for the root)
     */
    double getBranchLength(NodeRef node) {
        NodeRef parent = tree.getParent(node);
        if (parent == null) {
            return 0.0;
        }
        return (tree.getNodeHeight(parent) - tree.getNodeHeight(node)) * branchRateModel.getBranchRate(tree, node);
    }

    /**
     * @return log(expm1(x) / x), i.e. log of (1/R) int_0^R exp(x s / R) ds, stable for all x
     */
    static double logExpm1OverX(double x) {
        if (Math.abs(x) < 1E-8) {
            return 0.5 * x;
        }
        if (x > 30.0) {
            return x - Math.log(x) + Math.log1p(-Math.exp(-x));
        }
        return Math.log(Math.expm1(x) / x);
    }

    @Override
    public double[] getSiteLogLikelihoods() {
        if (!siteLogLikelihoodsKnown) {
            updateMAT();
            computeSiteLogLikelihoods();
        }
        return siteLogLikelihoods.clone();
    }

    @Override
    public void makeDirty() {
        allNodesDirty = true;
        treeChangedSinceCheck = true;
    }

    @Override
    public void storeState() {
        for (int i = 0; i < recomputedCount; i++) {
            recomputedSinceStore[recomputedNodes[i]] = false;
        }
        recomputedCount = 0;
        mutationStore.commit();
        storedLogLikelihood = logLikelihood;
        modelChangedSinceStore = false;
    }

    @Override
    public void restoreState() {
        // the terms recomputed since the store may have used the rejected branch lengths
        for (int i = 0; i < recomputedCount; i++) {
            int n = recomputedNodes[i];
            recomputedSinceStore[n] = false;
            markDirty(n);
        }
        recomputedCount = 0;
        mutationStore.revert();
        logLikelihood = storedLogLikelihood;
        siteLogLikelihoodsKnown = false;
        if (modelChangedSinceStore) {
            // cached terms may have been computed with the rejected model values
            modelVersion++;
        }
    }

    @Override
    protected void acceptState() {
        // nothing to do: the journal of changes is cleared by the next storeState()
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == tree) {
            // The tree likelihood hears of tree changes itself: only note which branch lengths may have
            // changed. Other objects (the height parameters, re-fired by the tree) come with a TreeChangedEvent.
            if (object instanceof TreeChangedEvent) {
                TreeChangedEvent event = (TreeChangedEvent) object;
                if (event.getParameter() == null || !event.isHeightChanged()) {
                    // (a change of a node's height parameter cannot change the topology)
                    treeChangedSinceCheck = true;
                }
                if (event.isNodeChanged() && event.getNode() != null) {
                    NodeRef node = event.getNode();
                    markDirty(node.getNumber());
                    for (int i = 0; i < tree.getChildCount(node); i++) {
                        markDirty(tree.getChild(node, i).getNumber());
                    }
                } else {
                    allNodesDirty = true;
                }
            }
            return;
        }
        if (model == branchRateModel) {
            // branch rates only affect branch lengths
            allNodesDirty = true;
        } else {
            modelVersion++;
            modelChangedSinceStore = true;
        }
        // index -1: TreeDataLikelihood would otherwise read the index as a node number
        fireModelChanged(object, -1);
    }

    @Override
    @SuppressWarnings("rawtypes")
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        modelVersion++;
        modelChangedSinceStore = true;
        fireModelChanged(variable, -1);
    }

    @Override
    public int getTraitCount() {
        return 1;
    }

    @Override
    public int getTraitDim() {
        return patternCount;
    }

    @Override
    public int getPartitionCat() {
        return 1;
    }

    @Override
    public RateRescalingScheme getRateRescalingScheme() {
        return null;
    }

    @Override
    public void setCallback(TreeDataLikelihood treeDataLikelihood) {
    }

    @Override
    public int vectorizeNodeOperations(List<NodeOperation> nodeOperations, int[] operations) {
        throw new UnsupportedOperationException("MAT likelihood does not use node operations");
    }

    @Override
    public void setComputePostOrderStatisticsOnly(boolean computePostOrderStatisticsOnly) {
    }

    @Override
    public boolean providesPostOrderStatisticsOnly() {
        return false;
    }

    @Override
    public PreOrderSettings getPreOrderSettings() {
        return PreOrderSettings.getDefault();
    }

    @Override
    public boolean getPreferGPU() {
        return false;
    }

    @Override
    public boolean getUseAmbiguities() {
        return false;
    }

    @Override
    public PartialsRescalingScheme getRescalingScheme() {
        return PartialsRescalingScheme.NONE;
    }

    @Override
    public boolean getDelayRescalingUntilUnderflow() {
        return false;
    }

    @Override
    public TreeTraversal.TraversalType getOptimalTraversalType() {
        return TreeTraversal.TraversalType.POST_ORDER;
    }

    @Override
    public long getTotalCalculationCount() {
        return calculationCount;
    }

    @Override
    public String getReport() {
        return NAME + ": " + patternCount + " patterns, " + mutationStore.getMutationCount() + " mutations, " +
                mutationStore.getAmbiguityCount() + " ambiguities, rebuilt " + rebuildCount + " times";
    }
}
