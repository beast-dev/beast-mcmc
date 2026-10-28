/*
 * MutationHistoryResampler.java
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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.operators.SubtreeMoveHook;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import dr.evomodel.treedatalikelihood.MAT.BlockSampler.Block;
import dr.evomodel.treedatalikelihood.MAT.BlockSampler.Conditional;
import dr.evomodel.treedatalikelihood.MAT.BlockSampler.UniformConditional;

/**
 * Resamples the ancestral states of a MAT jointly with a subtree prune-and-regraft move, as a
 * {@link SubtreeMoveHook} attached to the operator.
 *
 * With mutation times integrated out (or, for explicit hits, drawn separately; see below), a mutational
 * history is the state of every pattern at every internal node. When the subtree below X ({@code node})
 * moves, with its parent P, from beside S ({@code sibling}) to the branch above D ({@code destination}),
 * the states of a block are drawn from their exact full conditional given all other states (a Gibbs block
 * update), by forward filtering and backward sampling. The block covers the same states for the move and
 * its reverse, so the proposal densities of the old states (under the old topology) and of the new ones
 * (under the new topology) give the Hastings ratio, whatever topology move it accompanies.
 *
 * There are two blocks:
 * <ul>
 * <li>the path block: P and the internal nodes on the path between S and D in the pruned tree, at every
 * variable site, which contains every node whose ambiguity can change;</li>
 * <li>the local block ({@link LocalMove}): P at every variable site, and at each site only the nodes whose
 * data there come only from X (whose presence the move changes), as in the subtree moves of Varilly et
 * al. (Nature, 2026). Its cost does not grow with the distance moved.</li>
 * </ul>
 * The local block is the default in XML wherever the model allows it ({@link #supportsLocalBlock}): on big-3
 * it gave two to three times the ESS per minute of the path block and of BEAGLE, with the same posterior
 * (optimization_notes/mat_local_spr/benchmark).
 *
 * Without per-site rates, explicit hits or unequal pattern weights, sites are handled sparsely (see
 * {@link SparsePass}): a site with no events around the block has the same state a throughout, and its
 * block states stay unchanged with a probability that depends only on a and its rate class (its site rate
 * category, with a discrete site rate model whose categories are sampled). So the sites that change are
 * drawn by thinning (stepping through the variable sites with geometric jumps at the highest probability
 * of a change, and keeping each site visited with its own probability over that), and the unchanged ones
 * contribute n_ca log P(unchanged | c, a) for the number n_ca of them in each rate class and state. Sites whose only events
 * around the block are ambiguities are handled in the same way, grouped by their state and by which of
 * the block's neighbours are ambiguous there (and by rate class): there are typically many such sites
 * (runs of N) but few groups. A move then costs O(events around the block + changes), not O(sites).
 *
 * Branch history models with explicit hits ({@link BranchHistoryModel#hasExplicitHistories()}) use the
 * path block: the ancestral states are drawn with the hits summed out, then each site's hits on every
 * branch touching the block are drawn given its end states, and the probabilities of the old and new hits
 * enter the Hastings ratio.
 */
public class MutationHistoryResampler implements SubtreeMoveHook {

    /** Pattern classes relative to a block's neighbourhood (see {@link #getPatternClasses}) */
    public static final byte HAS_SUBSTITUTIONS = 0;
    public static final byte NO_EVENTS = 1;
    public static final byte AMBIGUITIES_ONLY = 2;
    public static final byte CONSTANT = 3;

    private final MATDataLikelihoodDelegate delegate;
    private final DefaultMutationStore store;
    private final BranchHistoryModel historyModel;
    private final Tree tree;
    private final int stateCount;
    private final int patternCount;
    private final double[] patternWeights;
    private final boolean localBlock;
    // the local block alone barely mixes the history away from the moves: it needs a MutationHistoryGibbsOperator
    private boolean hasGibbsOperator = false;
    private boolean gibbsWarningGiven = false;

    // the sites' rate classes: their site rate categories (0 with a single category), or none (-1) with
    // per-site rates; sparse passes count and group sites by (rate class, state)
    private final int rateClassCount;
    private final int[] variableSites;
    // scores and draws the states of a block, site by site
    private final BlockSampler blocks;
    // sites are handled sparsely (see SparsePass) without per-site rates, explicit hits or unequal weights
    private final boolean sparse;

    // the current path move (between beforeSubtreeMove and afterSubtreeMove)
    private int[] block = null;
    private int[] movedNodes = null;   // the nodes whose parents a subtree move changes
    private int[] candidateSites = null;
    private PathPass oldPathPass = null;
    private double logProposalOfOldHistories;
    // the current local move
    private LocalMove localMove = null;
    private int lastChangedUniformCount = 0;

    // per-site work arrays
    private final byte[] siteClasses;
    private final int[] siteStamps;    // candidate sites of the current path move are stamped with stamp
    private int stamp = 0;
    private final double[] oldSiteLogProbabilities;
    private final IndexSet changedSites;     // the sites at which the current move changed any state
    private final IndexSet passChangedSites; // those changed by the current sparse pass (see SparsePass)
    private final long[] siteBits;     // work bits of sortSites
    private final int[] oldMarks;      // the sites of a local move's old and new passes, and its hot sites,
    private final int[] newMarks;      // are marked with localStamp
    private final int[] hotMarks;
    private int localStamp = 0;

    // per-node work arrays and caches
    private final int[] nodeMarks;     // a node is marked if its entry equals nodeStamp
    private int nodeStamp = 0;
    private final int[][] ambiguousRuns;
    private final long[] ambiguousRunsVersion;
    // the summed weights of the constant patterns in each state (which every node has), and the weight of
    // every variable pattern (when sites are handled sparsely)
    private final double[] constantComposition;
    private final double variableWeight;
    private final long[] constantBits; // the constant patterns, as a bit set

    /**
     * A resampler with the path block (the XML parser uses the local block by default where it can).
     */
    public MutationHistoryResampler(MATDataLikelihoodDelegate delegate) {
        this(delegate, false);
    }

    /**
     * @param localBlock resample the local block (the regrafted parent, and the nodes whose data at a site
     *                   come only from the moved subtree) instead of the path block
     */
    public MutationHistoryResampler(MATDataLikelihoodDelegate delegate, boolean localBlock) {
        this.delegate = delegate;
        this.store = delegate.getMutableMutationStore();
        this.historyModel = delegate.getHistoryModel();
        this.tree = delegate.getTree();
        this.stateCount = delegate.getPatternList().getDataType().getStateCount();
        this.patternCount = delegate.getPatternList().getPatternCount();
        final int nodeCount = tree.getNodeCount();
        this.nodeMarks = new int[nodeCount];
        this.ambiguousRuns = new int[nodeCount][];
        this.ambiguousRunsVersion = new long[nodeCount];
        this.siteClasses = new byte[patternCount];
        this.siteStamps = new int[patternCount];
        this.oldSiteLogProbabilities = new double[patternCount];
        this.changedSites = new IndexSet(patternCount);
        this.passChangedSites = new IndexSet(patternCount);
        this.siteBits = new long[(patternCount + 63) / 64];
        this.oldMarks = new int[patternCount];
        this.newMarks = new int[patternCount];
        this.hotMarks = new int[patternCount];

        this.patternWeights = delegate.getPatternWeights();
        final double[] weights = patternWeights;
        // patterns grouped by weight, so that the block sampler can look conditionals up by weight class
        Map<Double, Integer> classes = new HashMap<Double, Integer>();
        final int[] weightClass = new int[patternCount];
        for (int p = 0; p < patternCount; p++) {
            Integer c = classes.get(weights[p]);
            if (c == null) {
                c = classes.size();
                classes.put(weights[p], c);
            }
            weightClass[p] = c;
        }
        final double[] classWeights = new double[classes.size()];
        for (Map.Entry<Double, Integer> e : classes.entrySet()) {
            classWeights[e.getValue()] = e.getKey();
        }

        int variableCount = 0;
        int[] variable = new int[patternCount];
        for (int p = 0; p < patternCount; p++) {
            if (!delegate.isConstantPattern(p)) {
                variable[variableCount++] = p;
            }
        }
        this.variableSites = Arrays.copyOf(variable, variableCount);
        this.sparse = supportsLocalBlock(delegate);
        this.constantComposition = new double[stateCount];
        for (int p = 0; p < patternCount; p++) {
            if (delegate.isConstantPattern(p)) {
                constantComposition[delegate.getPatternList().getPatternState(0, p)] += weights[p];
            }
        }
        this.variableWeight = variableCount > 0 ? weights[variableSites[0]] : 1.0;
        this.rateClassCount = delegate.getCategoryCount();
        this.blocks = new BlockSampler(delegate, weightClass, classWeights,
                variableCount > 0 ? weightClass[variableSites[0]] : 0, rateClassCount);
        this.constantBits = new long[(patternCount + 63) >>> 6];
        for (int p = 0; p < patternCount; p++) {
            if (delegate.isConstantPattern(p)) {
                constantBits[p >>> 6] |= 1L << p;
            }
        }
        this.localBlock = localBlock;
        if (localBlock && !sparse) {
            throw new IllegalArgumentException("The local block needs a branch history model without explicit " +
                    "hits, no per-site rate parameter (site rate categories are supported) and equal pattern weights");
        }
    }

    /**
     * @return true if sites can be handled sparsely, which the local block needs: a branch history model
     * without explicit hits, no per-site rates, and the same weight for every variable pattern
     */
    public static boolean supportsLocalBlock(MATDataLikelihoodDelegate delegate) {
        if (delegate.getHistoryModel().hasExplicitHistories() || delegate.hasSiteRates()) {
            return false;
        }
        final double[] weights = delegate.getPatternWeights();
        double weight = Double.NaN;
        for (int p = 0; p < delegate.getPatternList().getPatternCount(); p++) {
            if (!delegate.isConstantPattern(p)) {
                if (weight == weight && weights[p] != weight) {
                    return false;
                }
                weight = weights[p];
            }
        }
        // (with no variable sites there is nothing to resample)
        return weight == weight;
    }

    public MATDataLikelihoodDelegate getDelegate() {
        return delegate;
    }

    /**
     * Called by a {@link MutationHistoryGibbsOperator} that updates this resampler's histories.
     */
    void setGibbsOperator() {
        hasGibbsOperator = true;
    }

    @Override
    public void beforeSubtreeMove(NodeRef node, NodeRef sibling, NodeRef destination) {
        if (localBlock) {
            if (!hasGibbsOperator && !gibbsWarningGiven) {
                // the local block resamples a node's states only when a move regrafts next to it, so excess
                // mutations elsewhere are almost never removed and the tree parameters follow them
                // (at the info level, as BEAST prints warnings only when verbose)
                Logger.getLogger("dr.evomodel").info("WARNING: the mutation history resampler of " + tree.getId() +
                        " uses the local block, but no mutationHistoryGibbsOperator updates its histories. " +
                        "Without one the histories mix far too slowly for a correct posterior: add a " +
                        "<mutationHistoryGibbsOperator> with this resampler to the " +
                        "operators (with a weight of about 10% of the subtree moves' weight).");
                gibbsWarningGiven = true;
            }
            beforeLocalMove(node, sibling, destination);
        } else {
            beforePathMove(node, sibling, destination);
        }
    }

    @Override
    public double afterSubtreeMove() {
        return localBlock ? afterLocalMove() : afterPathMove();
    }

    /**
     * @return the number of sites with no events around the block that the last move changed (for testing)
     */
    public int getLastChangedUniformSiteCount() {
        return lastChangedUniformCount;
    }

    // ---------------------------------------------------------------------------------------------
    // The path block
    // ---------------------------------------------------------------------------------------------

    private void beforePathMove(NodeRef node, NodeRef sibling, NodeRef destination) {
        delegate.updateMAT();
        delegate.ensureRates();
        movedNodes = new int[]{sibling.getNumber(), tree.getParent(node).getNumber(), destination.getNumber()};
        block = getBlock(tree.getParent(node), sibling, destination);
        if (sparse) {
            oldPathPass = makePathPass(block, node.getNumber());
            sampleBlock(block, false, siteClasses, oldPathPass.substitutionSites, oldSiteLogProbabilities, oldPathPass);
        } else {
            candidateSites = getVariableSites(block);
            oldPathPass = null;
            sampleBlock(block, false, siteClasses, candidateSites, oldSiteLogProbabilities, null);
        }
        if (historyModel.hasExplicitHistories()) {
            // every candidate is resampled, so score all their hits around the block
            logProposalOfOldHistories = scoreHistories(getIncidentBranches(block), candidateSites);
        }
    }

    /**
     * Finds the sites with events around the block (setting candidateSites; the others are handled together,
     * see SparsePass) and makes the sparse pass for them.
     *
     * @param referenceNode a neighbour of the block whose states do not change (the uniform sites' states are
     *                      read there)
     */
    private PathPass makePathPass(int[] blockNodes, int referenceNode) {
        candidateSites = getCandidateSites(blockNodes);
        final int[] uniformCounts = getStateCounts(referenceNode);
        int substitutionCount = 0;
        int ambiguityCount = 0;
        int[] substitution = new int[candidateSites.length];
        int[] ambiguity = new int[candidateSites.length];
        for (int p : candidateSites) {
            if (delegate.isConstantPattern(p)) {
                continue;
            }
            int a = store.getState(referenceNode, p);
            if (a >= 0) {
                uniformCounts[countIndex(p, a)]--;
            }
            if (siteClasses[p] == AMBIGUITIES_ONLY) {
                ambiguity[ambiguityCount++] = p;
            } else {
                substitution[substitutionCount++] = p;
            }
        }
        // the ambiguity-only sites in pattern order, so that they can be swept against the neighbours' runs
        return new PathPass(referenceNode, uniformCounts, Arrays.copyOf(substitution, substitutionCount),
                sortSites(ambiguity, ambiguityCount), stamp);
    }

    // ---------------------------------------------------------------------------------------------
    // Gibbs updates in place
    // ---------------------------------------------------------------------------------------------

    /**
     * Draws the states of a connected block of internal nodes from their full conditional given all the other
     * states (and, with explicit hits, the hits on the branches touching the block given their end states),
     * without changing the topology: a Gibbs update of the mutational history (see
     * {@link MutationHistoryGibbsOperator}). The draw is that of the path block's new states.
     */
    public void resampleInPlace(int[] blockNodes) {
        delegate.updateMAT();
        delegate.ensureRates();

        final int[] resampled;
        if (sparse) {
            final PathPass pass = makePathPass(blockNodes, getOutsideChild(blockNodes));
            sampleBlock(blockNodes, true, siteClasses, pass.substitutionSites, null, pass);
            resampled = pass.substitutionSites;
            changedSites.addAll(pass.otherChangedSites);
        } else {
            resampled = getVariableSites(blockNodes);
            sampleBlock(blockNodes, true, siteClasses, resampled, null, null);
        }
        changedSites.addAll(resampled, resampled.length);
        delegate.setHistoryResampled();
        store.setFitchSetsStale(true);

        if (historyModel.hasExplicitHistories()) {
            // draw the hits of the resampled sites on every branch touching the block, given their end states
            List<Integer> branches = new ArrayList<Integer>();
            for (int n : getIncidentBranches(blockNodes)) {
                branches.add(n);
            }
            sampleHistories(branches, sortSites(changedSites));
            changedSites.clear();
        } else {
            // (no parent has changed)
            updateBranchesAroundBlock(blockNodes, null);
        }
        candidateSites = null;
        delegate.fireHistoryChanged();
    }

    /**
     * @return a child of one of the block's nodes that is not in the block
     */
    private int getOutsideChild(int[] blockNodes) {
        for (int n : blockNodes) {
            NodeRef node = tree.getNode(n);
            for (int i = 0; i < tree.getChildCount(node); i++) {
                int c = tree.getChild(node, i).getNumber();
                if (!contains(blockNodes, c)) {
                    return c;
                }
            }
        }
        throw new IllegalArgumentException("The block has no children outside it");
    }

    private double afterPathMove() {
        if (block == null) {
            throw new IllegalStateException("afterSubtreeMove called without beforeSubtreeMove");
        }
        final int[] currentBlock = block;
        final PathPass oldPass = oldPathPass;
        block = null;
        oldPathPass = null;
        final IntBuffer reparented = syncParents(movedNodes);

        // every candidate is resampled (the others are handled together when sampling sparsely)
        final int[] resampled = oldPass != null ? oldPass.substitutionSites : candidateSites;
        changedSites.addAll(resampled, resampled.length);

        double logProposalOfOldStates = oldPass != null ?
                oldPass.uniformLogProbability + oldPass.ambiguityLogProbability : 0.0;
        for (int p : resampled) {
            logProposalOfOldStates += oldSiteLogProbabilities[p];
        }
        final PathPass newPass = oldPass == null ? null :
                new PathPass(oldPass.referenceNode, oldPass.uniformCounts, oldPass.substitutionSites,
                        oldPass.ambiguitySites, oldPass.candidateStamp);
        double logProposalOfNewStates = sampleBlock(currentBlock, true, siteClasses, resampled, null, newPass);
        lastChangedUniformCount = newPass == null ? 0 : newPass.changedUniformCount;
        delegate.setHistoryResampled();

        if (newPass != null) {
            // the other sites whose block states changed also need their branches updated
            changedSites.addAll(newPass.otherChangedSites);
        }

        // the Fitch sets no longer match the topology; parsimony will recompute them if it is ever needed
        store.setFitchSetsStale(true);

        if (historyModel.hasExplicitHistories()) {
            // draw the hits of the resampled sites on every branch touching the block, given their end states
            Set<Integer> branches = new HashSet<Integer>();
            for (int i = 0; i < reparented.count; i++) {
                branches.add(reparented.values[i]);
            }
            for (int n : getIncidentBranches(currentBlock)) {
                branches.add(n);
            }
            double logProposalOfNewHistories = sampleHistories(branches, resampled);
            changedSites.clear();
            candidateSites = null;
            delegate.markTopologySynchronised();
            return logProposalOfOldStates - logProposalOfNewStates + logProposalOfOldHistories - logProposalOfNewHistories;
        }

        updateBranchesAroundBlock(currentBlock, reparented);
        candidateSites = null;
        delegate.markTopologySynchronised();

        return logProposalOfOldStates - logProposalOfNewStates;
    }

    /**
     * Recomputes the events of the branches around a block whose states have changed at changedSites (which
     * is then emptied): in full for the nodes whose parent changed, and at the changed sites only for the
     * block's other nodes and their children.
     *
     * @param reparented the nodes whose parent changed, or null
     */
    private void updateBranchesAroundBlock(int[] blockNodes, IntBuffer reparented) {
        final int[] sortedChanges = sortSites(changedSites);
        final int touchedStamp = nextNodeStamp();
        if (reparented != null) {
            for (int i = 0; i < reparented.count; i++) {
                int n = reparented.values[i];
                nodeMarks[n] = touchedStamp;
                store.updateBranch(n);
            }
        }
        final IntBuffer touched = new IntBuffer();
        for (int n : blockNodes) {
            addWithChildren(touched, n, touchedStamp);
        }
        for (int i = 0; i < touched.count; i++) {
            store.updateBranchSites(touched.values[i], changedSites, sortedChanges);
        }
        changedSites.clear();
    }

    /**
     * @return the regrafted parent plus the internal nodes on the path between a and b, in the order of the
     * set (which fixes the order in which the candidate sites are found, and so drawn)
     */
    public int[] getBlock(NodeRef parent, NodeRef a, NodeRef b) {
        Set<NodeRef> ancestorsOfA = new HashSet<NodeRef>();
        for (NodeRef n = a; n != null; n = tree.getParent(n)) {
            ancestorsOfA.add(n);
        }
        List<NodeRef> path = new ArrayList<NodeRef>();
        NodeRef mrca = b;
        while (!ancestorsOfA.contains(mrca)) {
            path.add(mrca);
            mrca = tree.getParent(mrca);
        }
        for (NodeRef n = a; n != mrca; n = tree.getParent(n)) {
            path.add(n);
        }
        path.add(mrca);

        Set<Integer> nodes = new HashSet<Integer>();
        nodes.add(parent.getNumber());
        for (NodeRef n : path) {
            if (!tree.isExternal(n)) {
                nodes.add(n.getNumber());
            }
        }
        int[] result = new int[nodes.size()];
        int i = 0;
        for (int n : nodes) {
            result[i++] = n;
        }
        return result;
    }

    /**
     * @return the sites with an event (substitution or ambiguity) on a branch within the block's
     * neighbourhood (the branches above the block nodes and their children), with their classes set in
     * siteClasses and stamped in siteStamps. Every other site has the same state throughout the neighbourhood.
     */
    int[] getCandidateSites(int[] blockNodes) {
        stamp++;
        int[] sites = new int[16];
        int count = 0;
        for (int n : blockNodes) {
            NodeRef node = tree.getNode(n);
            for (int i = -1; i < tree.getChildCount(node); i++) {
                BranchMutations mutations = store.getBranchMutations(i < 0 ? n : tree.getChild(node, i).getNumber());
                for (int k = 0; k < mutations.size(); k++) {
                    int p = mutations.getPattern(k);
                    if (siteStamps[p] != stamp) {
                        siteStamps[p] = stamp;
                        siteClasses[p] = mutations.isAmbiguity(k) ? AMBIGUITIES_ONLY : HAS_SUBSTITUTIONS;
                        if (count == sites.length) {
                            sites = Arrays.copyOf(sites, 2 * count);
                        }
                        sites[count++] = p;
                    } else if (!mutations.isAmbiguity(k)) {
                        siteClasses[p] = HAS_SUBSTITUTIONS;
                    }
                }
            }
        }
        return Arrays.copyOf(sites, count);
    }

    /**
     * @return every non-constant site, with its class around the block set in siteClasses
     */
    private int[] getVariableSites(int[] blockNodes) {
        byte[] classes = getPatternClasses(blockNodes);
        int[] sites = new int[patternCount];
        int count = 0;
        for (int p = 0; p < patternCount; p++) {
            if (!delegate.isConstantPattern(p)) {
                siteClasses[p] = classes[p];
                sites[count++] = p;
            }
        }
        return Arrays.copyOf(sites, count);
    }

    /**
     * Classifies the patterns by the events on the branches within the block's neighbourhood (the
     * block, its children and the top node's parent). With no substitutions there, every node in the
     * neighbourhood that is not ambiguous has the same state. The neighbourhood after the move is
     * contained in the neighbourhood before it, so this also holds in the new topology.
     *
     * Constant patterns are assumed to have no mutations and are classed CONSTANT.
     *
     * @return for each pattern HAS_SUBSTITUTIONS, NO_EVENTS, AMBIGUITIES_ONLY or CONSTANT
     */
    public byte[] getPatternClasses(int[] blockNodes) {
        byte[] classes = new byte[patternCount];
        for (int p = 0; p < patternCount; p++) {
            classes[p] = delegate.isConstantPattern(p) ? CONSTANT : NO_EVENTS;
        }
        for (int n : blockNodes) {
            classify(store.getBranchMutations(n), classes);
            NodeRef node = tree.getNode(n);
            for (int i = 0; i < tree.getChildCount(node); i++) {
                classify(store.getBranchMutations(tree.getChild(node, i).getNumber()), classes);
            }
        }
        return classes;
    }

    private static void classify(BranchMutations mutations, byte[] classes) {
        for (int k = 0; k < mutations.size(); k++) {
            int p = mutations.getPattern(k);
            if (!mutations.isAmbiguity(k)) {
                classes[p] = HAS_SUBSTITUTIONS;
            } else if (classes[p] == NO_EVENTS) {
                classes[p] = AMBIGUITIES_ONLY;
            }
        }
    }

    private static boolean contains(int[] nodes, int n) {
        for (int b : nodes) {
            if (b == n) {
                return true;
            }
        }
        return false;
    }

    // ---------------------------------------------------------------------------------------------
    // The local block
    // ---------------------------------------------------------------------------------------------

    /**
     * What a local move carries from before the topology change to after it.
     * <p>
     * The subtree below X ({@code node}) moves, with its parent P, from beside S ({@code sibling}) to the
     * branch above D ({@code destination}). Before the move P's neighbours are its parent G, S and X; after
     * it they are its new parent N, D and X. P's states are resampled at every variable site, given its
     * neighbours. At a site where some of X's ancestors have data only from X (a hot site), their states
     * change too: on the old side G, and each ancestor above it whose other children have no data there,
     * lose their states (which is determined by the topology); on the new side N, and the ancestors above it
     * that had no data there, gain states, which are drawn with P's as a chain from the first ancestor that
     * had data down to X. No other state changes. A move and its reverse change the same states and each
     * draws them from their full conditional given all the others, so the Hastings ratio is the ratio of
     * the old states' conditional (under the old topology) to the new states' (under the new one).
     * <p>
     * The sites are found from the substitutions on the few branches that relate P's neighbours and from the
     * runs of ambiguous patterns of the neighbours (see {@link LocalPass}), so a move costs
     * O(events near P + runs + changes), whatever the distance moved
     * (see optimization_notes/mat_local_spr/DESIGN.md).
     */
    private static final class LocalMove {
        int node;                      // X
        int parent;                    // P
        int sibling;                   // S
        int destination;               // D
        int newParent;                 // N: -1 if P becomes the root
        int[] presentX;                // runs of patterns at which X has data
        int[] oldChain;                // G and the ancestors above it, as far as any site's old chain reaches
        int[][] oldChainLevels;        // [i]: runs of patterns at which the old chain reaches oldChain[i]
        int[] oldHotSites;             // the variable sites at which the old chain is not empty
        int[] newSubstitutionSites;
        int newStamp;
        double logProposalOfOldStates;
    }

    private void beforeLocalMove(NodeRef node, NodeRef sibling, NodeRef destination) {
        delegate.updateMAT();
        delegate.ensureRates();
        localStamp++;
        final int markStamp = localStamp;

        final LocalMove move = new LocalMove();
        final NodeRef pNode = tree.getParent(node);
        final NodeRef gNode = tree.getParent(pNode);
        final NodeRef dParent = tree.getParent(destination);
        // P's new parent: D's parent, unless D is S (which then stays below G)
        final NodeRef newParentNode = dParent != null && dParent.getNumber() == pNode.getNumber() ? gNode : dParent;
        move.node = node.getNumber();
        move.parent = pNode.getNumber();
        move.sibling = sibling.getNumber();
        move.destination = destination.getNumber();
        move.newParent = newParentNode == null ? -1 : newParentNode.getNumber();
        move.presentX = SiteRuns.complement(getAmbiguousRuns(move.node), patternCount);
        findOldChain(move, gNode, newParentNode);
        for (int site : move.oldHotSites) {
            hotMarks[site] = markStamp;
        }

        // The old pass: P with neighbours G, S and X.
        // (the sites at which a neighbour is ambiguous are found from their runs when P is scored or drawn)
        final IntBuffer oldSubstitutions = new IntBuffer();
        collectSubstitutionSites(gNode == null ? new int[]{move.node, move.sibling} :
                new int[]{move.node, move.sibling, move.parent}, oldMarks, markStamp, oldSubstitutions);
        final int[] oldSubstitutionSites = oldSubstitutions.toArray();
        final int[] oldUniformCounts = getUniformCounts(move.node, oldSubstitutionSites, move.oldHotSites, oldMarks, markStamp);

        // The new pass: P with neighbours N, D and X. P's state there is affected by the substitutions
        // between X and N (in the old topology, where they are on the path between them) and between D and
        // N, and by the runs of ambiguous patterns of X, D and N (after the old chain loses its states, since
        // D may be in it).
        final NodeRef anchor = newParentNode != null ? newParentNode : destination;
        final IntBuffer pathBranches = getBranchesBetween(node, anchor);
        if (newParentNode != null) {
            pathBranches.add(move.destination);
        }
        final IntBuffer newSubstitutions = new IntBuffer();
        collectSubstitutionSites(pathBranches.toArray(), newMarks, markStamp, newSubstitutions);
        move.newSubstitutionSites = newSubstitutions.toArray();
        move.newStamp = markStamp;

        // score P's current states (and the old chain's) under the old topology
        final LocalPass pass = new LocalPass(move.node, oldUniformCounts, oldSubstitutionSites, move.oldHotSites);
        final int[] substitutionSites = withoutMarked(oldSubstitutionSites, hotMarks, markStamp);
        double logProbability = sampleBlock(new int[]{move.parent}, false, null, substitutionSites, null, pass);
        logProbability += sampleChains(move.parent, move.oldChain, move.oldChainLevels, false);
        move.logProposalOfOldStates = logProbability;
        localMove = move;
    }

    /**
     * Finds the old chain: G and each ancestor above it lose their states where X was their only data, up to
     * the first node that is above X's new position too.
     */
    private void findOldChain(LocalMove move, NodeRef gNode, NodeRef newParentNode) {
        final int newAncestorStamp = nextNodeStamp();
        for (NodeRef n = newParentNode; n != null; n = tree.getParent(n)) {
            nodeMarks[n.getNumber()] = newAncestorStamp;
        }
        final IntBuffer chain = new IntBuffer();
        final List<int[]> levels = new ArrayList<int[]>();
        int[] hot = SiteRuns.intersect(move.presentX, getAmbiguousRuns(move.sibling));
        int previous = move.parent;
        for (NodeRef n = gNode; n != null && hot.length > 0 && nodeMarks[n.getNumber()] != newAncestorStamp; n = tree.getParent(n)) {
            for (int i = 0; i < tree.getChildCount(n); i++) {
                int c = tree.getChild(n, i).getNumber();
                if (c != previous) {
                    hot = SiteRuns.intersect(hot, getAmbiguousRuns(c));
                }
            }
            if (hot.length == 0) {
                break;
            }
            chain.add(n.getNumber());
            levels.add(hot);
            previous = n.getNumber();
        }
        move.oldChain = chain.toArray();
        move.oldChainLevels = levels.toArray(new int[levels.size()][]);
        move.oldHotSites = levels.isEmpty() ? new int[0] : variableSitesIn(levels.get(0));
    }

    /**
     * @return the nodes whose branches lie on the path between two nodes: those from b up to (not including)
     * the first ancestor of a, then those from a up to it
     */
    private IntBuffer getBranchesBetween(NodeRef a, NodeRef b) {
        final int ancestorStamp = nextNodeStamp();
        for (NodeRef n = a; n != null; n = tree.getParent(n)) {
            nodeMarks[n.getNumber()] = ancestorStamp;
        }
        final IntBuffer branches = new IntBuffer();
        NodeRef mrca = b;
        while (nodeMarks[mrca.getNumber()] != ancestorStamp) {
            branches.add(mrca.getNumber());
            mrca = tree.getParent(mrca);
        }
        for (NodeRef n = a; n.getNumber() != mrca.getNumber(); n = tree.getParent(n)) {
            branches.add(n.getNumber());
        }
        return branches;
    }

    private double afterLocalMove() {
        final LocalMove move = localMove;
        if (move == null) {
            throw new IllegalStateException("afterSubtreeMove called without beforeSubtreeMove");
        }
        localMove = null;
        final IntBuffer reparented = syncParents(new int[]{move.sibling, move.parent, move.destination});

        // The new chain: N, and each ancestor above it that had no data, gain states where X has data
        // (read before the old chain loses its states, so as states before the move).
        final IntBuffer chain = new IntBuffer();
        final List<int[]> levels = new ArrayList<int[]>();
        if (move.newParent >= 0) {
            int[] hot = SiteRuns.intersect(move.presentX, getAmbiguousRuns(move.newParent));
            NodeRef n = tree.getNode(move.newParent);
            while (hot.length > 0) {
                chain.add(n.getNumber());
                levels.add(hot);
                NodeRef parent = tree.getParent(n);
                if (parent == null) {
                    break;
                }
                hot = SiteRuns.intersect(hot, getAmbiguousRuns(parent.getNumber()));
                n = parent;
            }
        }
        final int[] newChain = chain.toArray();
        final int[][] newChainLevels = levels.toArray(new int[levels.size()][]);
        final int[] newHotSites = levels.isEmpty() ? new int[0] : variableSitesIn(levels.get(0));
        localStamp++;
        final int newHotStamp = localStamp;
        for (int site : newHotSites) {
            hotMarks[site] = newHotStamp;
        }

        // the old chain loses its states at its hot sites
        for (int i = 0; i < move.oldChain.length; i++) {
            for (int site : variableSitesIn(move.oldChainLevels[i])) {
                store.setState(move.oldChain[i], site, -1);
            }
        }

        // draw P's states (and the new chain's) under the new topology
        final LocalPass pass = new LocalPass(move.node,
                getUniformCounts(move.node, move.newSubstitutionSites, newHotSites, newMarks, move.newStamp),
                move.newSubstitutionSites, newHotSites);
        final int[] substitutionSites = withoutMarked(move.newSubstitutionSites, hotMarks, newHotStamp);
        double logProposalOfNewStates = sampleBlock(new int[]{move.parent}, true, null, substitutionSites, null, pass);
        logProposalOfNewStates += sampleChains(move.parent, newChain, newChainLevels, true);
        lastChangedUniformCount = pass.changedUniformCount;
        delegate.setHistoryResampled();
        // the Fitch sets no longer match the topology; parsimony will recompute them if it is ever needed
        store.setFitchSetsStale(true);

        // the sites at which any state changed
        changedSites.addAll(substitutionSites, substitutionSites.length);
        changedSites.addAll(pass.otherChangedSites);
        changedSites.addAll(move.oldHotSites, move.oldHotSites.length);
        changedSites.addAll(newHotSites, newHotSites.length);

        // the nodes whose states changed: P (whose children now include X) and the chains
        final int[] changedNodes = new int[1 + move.oldChain.length + newChain.length];
        changedNodes[0] = move.parent;
        System.arraycopy(move.oldChain, 0, changedNodes, 1, move.oldChain.length);
        System.arraycopy(newChain, 0, changedNodes, 1 + move.oldChain.length, newChain.length);
        updateBranchesAroundBlock(changedNodes, reparented);
        delegate.markTopologySynchronised();

        return move.logProposalOfOldStates - logProposalOfNewStates;
    }

    /**
     * Scores or draws the states of the chains {P, chain[0], ..., chain[k]} at the sites where the chain
     * has length k + 1 (levels[i] are the runs of patterns at which it reaches chain[i]).
     */
    private double sampleChains(int p, int[] chain, int[][] levels, boolean sample) {
        double logProbability = 0.0;
        for (int k = 0; k < chain.length; k++) {
            int[] exact = k + 1 < chain.length ? SiteRuns.subtract(levels[k], levels[k + 1], patternCount) : levels[k];
            int[] sites = variableSitesIn(exact);
            if (sites.length == 0) {
                continue;
            }
            int[] chainBlock = new int[k + 2];
            chainBlock[0] = p;
            System.arraycopy(chain, 0, chainBlock, 1, k + 1);
            logProbability += sampleBlock(chainBlock, sample, null, sites, null, null);
        }
        return logProbability;
    }

    /**
     * Adds the variable sites with a substitution on the branches above the given nodes (and not yet marked).
     */
    private void collectSubstitutionSites(int[] branches, int[] marks, int markStamp, IntBuffer out) {
        for (int n : branches) {
            BranchMutations mutations = store.getBranchMutations(n);
            for (int k = 0; k < mutations.getSubstitutionCount(); k++) {
                int site = mutations.getPattern(k);
                if (marks[site] != markStamp && !delegate.isConstantPattern(site)) {
                    marks[site] = markStamp;
                    out.add(site);
                }
            }
        }
    }

    /**
     * @return the number of variable sites in each state at node x, excluding the substitution sites (marked)
     * and the hot sites (the sites at which a neighbour of P is ambiguous are excluded when they are swept)
     */
    private int[] getUniformCounts(int x, int[] substitutions, int[] hot, int[] marks, int markStamp) {
        final int[] counts = getStateCounts(x);
        for (int site : substitutions) {
            int a = store.getState(x, site);
            if (a >= 0) {
                counts[countIndex(site, a)]--;
            }
        }
        for (int site : hot) {
            int a = store.getState(x, site);
            if (a >= 0 && marks[site] != markStamp) {
                counts[countIndex(site, a)]--;
            }
        }
        return counts;
    }

    private static int[] withoutMarked(int[] sites, int[] marks, int markStamp) {
        int[] result = new int[sites.length];
        int count = 0;
        for (int site : sites) {
            if (marks[site] != markStamp) {
                result[count++] = site;
            }
        }
        return count == sites.length ? sites : Arrays.copyOf(result, count);
    }

    /**
     * @return the variable sites in the runs, in pattern order
     */
    private int[] variableSitesIn(int[] runs) {
        IntBuffer sites = new IntBuffer();
        for (int i = 0; i < runs.length; i += 2) {
            for (int site = runs[i]; site < runs[i + 1]; site++) {
                if (!delegate.isConstantPattern(site)) {
                    sites.add(site);
                }
            }
        }
        return sites.toArray();
    }

    private void addWithChildren(IntBuffer nodes, int n, int markStamp) {
        addUnmarked(nodes, n, markStamp);
        NodeRef node = tree.getNode(n);
        for (int i = 0; i < tree.getChildCount(node); i++) {
            addUnmarked(nodes, tree.getChild(node, i).getNumber(), markStamp);
        }
    }

    private void addUnmarked(IntBuffer nodes, int n, int markStamp) {
        if (nodeMarks[n] != markStamp) {
            nodeMarks[n] = markStamp;
            nodes.add(n);
        }
    }

    private int nextNodeStamp() {
        return ++nodeStamp;
    }

    /**
     * Gives the store the tree's parents of the nodes a subtree move reparents (the only ones it changes).
     *
     * @return the nodes whose parent changed
     */
    private IntBuffer syncParents(int[] nodes) {
        final IntBuffer reparented = new IntBuffer();
        for (int n : nodes) {
            NodeRef parent = tree.getParent(tree.getNode(n));
            int parentNumber = parent == null ? -1 : parent.getNumber();
            if (store.getParent(n) != parentNumber) {
                store.setParent(n, parentNumber);
                reparented.add(n);
            }
        }
        if (store.getRoot() < 0 || store.getParent(store.getRoot()) >= 0) {
            throw new IllegalStateException("The subtree move changed more parents than expected");
        }
        return reparented;
    }

    // ---------------------------------------------------------------------------------------------
    // Whole-tree updates and explicit hits
    // ---------------------------------------------------------------------------------------------

    /**
     * @return the log probability of all the ancestral states (of the variable sites) and hits under their
     * full conditional distribution given the tree; the likelihood minus this is the likelihood with
     * them summed out
     */
    public double getLogFullConditional() {
        delegate.updateMAT();
        delegate.ensureRates();
        double logProbability = sampleBlock(getInternalNodes(), false);
        if (historyModel.hasExplicitHistories()) {
            logProbability += scoreHistories(getNonRootNodes(), getVariableSites(getInternalNodes()));
        }
        return logProbability;
    }

    /**
     * Draws all the ancestral states (of the variable sites) and hits from their full conditional
     * distribution given the tree (a Gibbs update of the whole MAT).
     */
    public void resampleAll() {
        delegate.updateMAT();
        delegate.ensureRates();
        int[] internal = getInternalNodes();
        sampleBlock(internal, true);
        if (historyModel.hasExplicitHistories()) {
            int[] sites = getVariableSites(internal);
            changedSites.addAll(sites, sites.length);
            List<Integer> branches = new ArrayList<Integer>();
            for (int n : getNonRootNodes()) {
                branches.add(n);
            }
            sampleHistories(branches, sites);
            changedSites.clear();
        } else {
            store.updateAllBranches();
        }
        delegate.setHistoryResampled();
    }

    private int[] getInternalNodes() {
        int[] nodes = new int[tree.getInternalNodeCount()];
        for (int i = 0; i < nodes.length; i++) {
            nodes[i] = tree.getInternalNode(i).getNumber();
        }
        return nodes;
    }

    private int[] getNonRootNodes() {
        int[] nodes = new int[tree.getNodeCount() - 1];
        int i = 0;
        for (int n = 0; n < tree.getNodeCount(); n++) {
            if (!tree.isRoot(tree.getNode(n))) {
                nodes[i++] = n;
            }
        }
        return nodes;
    }

    /**
     * @return the nodes whose branches touch the block in the current topology: the block nodes (other
     * than the root) and their children
     */
    private int[] getIncidentBranches(int[] blockNodes) {
        Set<Integer> branches = new HashSet<Integer>();
        for (int n : blockNodes) {
            NodeRef node = tree.getNode(n);
            if (!tree.isRoot(node)) {
                branches.add(n);
            }
            for (int i = 0; i < tree.getChildCount(node); i++) {
                branches.add(tree.getChild(node, i).getNumber());
            }
        }
        int[] result = new int[branches.size()];
        int i = 0;
        for (int n : branches) {
            result[i++] = n;
        }
        return result;
    }

    /**
     * @return the log probability of the current hits of the given sites on the branches above the given
     * nodes, given their end states
     */
    private double scoreHistories(int[] branches, int[] sites) {
        final double[] rates = delegate.getPatternRates();
        final int[] start = new int[patternCount];
        final int[] end = new int[patternCount];
        double logProbability = 0.0;
        for (int n : branches) {
            int parent = store.getParent(n);
            double length = delegate.getBranchLength(tree.getNode(n));
            BranchMutations mutations = store.getBranchMutations(n);
            for (int p : sites) {
                start[p] = -1;
            }
            int k = 0;
            while (k < mutations.getSubstitutionCount()) {
                int siteEnd = mutations.getSiteEnd(k);
                start[mutations.getPattern(k)] = k;
                end[mutations.getPattern(k)] = siteEnd;
                k = siteEnd;
            }
            for (int p : sites) {
                int b = store.getState(n, p);
                if (b < 0) {
                    continue;
                }
                int a = store.getState(parent, p);
                logProbability += start[p] < 0 ?
                        historyModel.getHistoryLogProposal(mutations, 0, 0, a, b, length, rates[p]) :
                        historyModel.getHistoryLogProposal(mutations, start[p], end[p], a, b, length, rates[p]);
            }
        }
        return logProbability;
    }

    /**
     * Draws the hits of the sites in changedSites (listed in pattern order in sites) on the branches above
     * the given nodes, given their end states, keeping the events of the other sites.
     *
     * @return the log probability of the drawn hits
     */
    private double sampleHistories(Iterable<Integer> branches, int[] sites) {
        final double[] rates = delegate.getPatternRates();
        double logProbability = 0.0;
        for (int n : branches) {
            int parent = store.getParent(n);
            if (parent < 0) {
                continue;
            }
            double length = delegate.getBranchLength(tree.getNode(n));
            BranchMutations old = store.getBranchMutations(n);
            BranchMutations.Builder builder = new BranchMutations.Builder(true);
            int k = 0;
            while (k < old.getSubstitutionCount()) {
                int siteEnd = old.getSiteEnd(k);
                if (!changedSites.contains(old.getPattern(k))) {
                    for (int i = k; i < siteEnd; i++) {
                        builder.addSubstitution(old.getPattern(i), old.getFromState(i), old.getToState(i),
                                old.hasTimes() ? old.getTime(i) : (i - k + 1.0) / (siteEnd - k + 1.0));
                    }
                }
                k = siteEnd;
            }
            for (int i = old.getSubstitutionCount(); i < old.size(); i++) {
                if (!changedSites.contains(old.getPattern(i))) {
                    builder.addAmbiguity(old.getPattern(i), old.getFromState(i));
                }
            }
            for (int p : sites) {
                int a = store.getState(parent, p);
                int b = store.getState(n, p);
                if (b < 0) {
                    if (a >= 0) {
                        builder.addAmbiguity(p, a);
                    }
                } else {
                    logProbability += historyModel.sampleHistory(p, a, b, length, rates[p], builder);
                }
            }
            store.setBranchMutations(n, builder.build());
        }
        return logProbability;
    }

    // ---------------------------------------------------------------------------------------------
    // Per-node caches and site lists
    // ---------------------------------------------------------------------------------------------

    /**
     * @return the runs of patterns at which the node is ambiguous (cached for the node's current states)
     */
    private int[] getAmbiguousRuns(int n) {
        long version = store.getStateVersion(n);
        if (ambiguousRuns[n] == null || ambiguousRunsVersion[n] != version) {
            ambiguousRuns[n] = SiteRuns.ambiguous(store.getStates(n), patternCount);
            ambiguousRunsVersion[n] = version;
        }
        return ambiguousRuns[n];
    }

    /**
     * @return the number of the variable sites in each state at the node (when sites are handled sparsely,
     * so that every variable pattern has the same weight), from the node's composition: constant patterns
     * have their state at every node
     */
    private int[] getStateCounts(int n) {
        final int[] counts = new int[rateClassCount * stateCount];
        if (rateClassCount > 1) {
            // [category * stateCount + state], from the node's category composition (variable patterns only)
            final double[] composition = store.getCategoryComposition(n);
            for (int i = 0; i < counts.length; i++) {
                counts[i] = (int) Math.round(composition[i] / variableWeight);
            }
        } else {
            final double[] composition = store.getComposition(n);
            for (int a = 0; a < stateCount; a++) {
                counts[a] = (int) Math.round((composition[a] - constantComposition[a]) / variableWeight);
            }
        }
        assert Arrays.equals(counts, countStates(n)) : "state counts of node " + n;
        return counts;
    }

    private int[] countStates(int n) {
        final int[] counts = new int[rateClassCount * stateCount];
        final byte[] states = store.getStates(n);
        for (int p : variableSites) {
            if (states[p] >= 0) {
                counts[countIndex(p, states[p])]++;
            }
        }
        return counts;
    }

    /**
     * @return the index of a site in state a in the counts of the sparse passes: (rate class, state)
     */
    private int countIndex(int p, int a) {
        return rateClassCount > 1 ? store.getCategory(p) * stateCount + a : a;
    }

    /**
     * @return the given sites in pattern order (each site at most once)
     */
    private int[] sortSites(int[] sites, int count) {
        for (int i = 0; i < count; i++) {
            siteBits[sites[i] >> 6] |= 1L << (sites[i] & 63);
        }
        return takeSortedSites(count);
    }

    private int[] sortSites(IndexSet sites) {
        for (int i = 0; i < sites.size(); i++) {
            siteBits[sites.get(i) >> 6] |= 1L << (sites.get(i) & 63);
        }
        return takeSortedSites(sites.size());
    }

    /**
     * @return the sites set in siteBits (which it clears), in pattern order
     */
    private int[] takeSortedSites(int count) {
        int[] sorted = new int[count];
        int k = 0;
        for (int word = 0; word < siteBits.length; word++) {
            long bits = siteBits[word];
            siteBits[word] = 0L;
            while (bits != 0L) {
                sorted[k++] = (word << 6) + Long.numberOfTrailingZeros(bits);
                bits &= bits - 1;
            }
        }
        return k == count ? sorted : Arrays.copyOf(sorted, k);
    }

    // ---------------------------------------------------------------------------------------------
    // Blocks (see BlockSampler)
    // ---------------------------------------------------------------------------------------------

    /**
     * Draws (or, if sample is false, scores) the states of the block's nodes from their full
     * conditional distribution given all other states in the current tree.
     *
     * @return the log probability of the drawn (or existing) states under that distribution
     */
    public double sampleBlock(int[] blockNodes, boolean sample) {
        return sampleBlock(blockNodes, sample, null, null, null, null);
    }

    /**
     * As {@link #sampleBlock(int[], boolean)} but with patterns that have no substitutions around the
     * block (see {@link #getPatternClasses}) handled by a fast exact path.
     */
    public double sampleBlock(int[] blockNodes, boolean sample, byte[] patternClasses) {
        return sampleBlock(blockNodes, sample, patternClasses, null, null, null);
    }

    /**
     * As {@link #sampleBlock(int[], boolean)} but only for the given sites (all sites if null), with the
     * patterns that have no substitutions around the block (see {@link #getPatternClasses}) handled by a
     * fast exact path, optionally recording each site's log probability, and then, if a sparse pass is
     * given, for the sites it handles (those with no events, or only ambiguities, around the block).
     */
    private double sampleBlock(int[] blockNodes, boolean sample, byte[] patternClasses,
                               int[] sites, double[] siteLogProbabilities, SparsePass pass) {
        delegate.ensureRates();
        final Block block = blocks.newBlock(blockNodes, sample);
        double logProbability = 0.0;
        final int siteCount = sites == null ? patternCount : sites.length;
        for (int siteIndex = 0; siteIndex < siteCount; siteIndex++) {
            final int p = sites == null ? siteIndex : sites[siteIndex];
            // constant sites are assumed to have no mutations: their states are never changed
            final double before = logProbability;
            if (!delegate.isConstantPattern(p)) {
                logProbability = blocks.sampleSite(block, p, patternClasses == null ? HAS_SUBSTITUTIONS : patternClasses[p],
                        logProbability);
            }
            if (siteLogProbabilities != null) {
                siteLogProbabilities[p] = logProbability - before;
            }
        }
        if (pass != null) {
            pass.run(block);
            logProbability += pass.ambiguityLogProbability;
            logProbability += pass.uniformLogProbability;
        }
        block.release();
        return logProbability;
    }

    // ---------------------------------------------------------------------------------------------
    // Sparse passes: the sites with no events, or only ambiguities, around the block
    // ---------------------------------------------------------------------------------------------

    /**
     * The sparse part of scoring or drawing a block (see the class comment): after the sites handled
     * individually, the sites whose only events around the block are ambiguities are handled in groups, then
     * the uniform sites (no events around the block) by their counts and thinning.
     */
    private abstract class SparsePass {
        final int referenceNode;   // X: the uniform sites' states are read here (the move does not change them)
        final int[] uniformCounts; // the number of uniform sites in each state
        double ambiguityLogProbability;
        double uniformLogProbability;
        // the sites, other than those handled individually, whose block states changed (when drawing)
        final IndexSet otherChangedSites = passChangedSites;
        int changedUniformCount;

        SparsePass(int referenceNode, int[] uniformCounts) {
            this.referenceNode = referenceNode;
            this.uniformCounts = uniformCounts;
        }

        /**
         * @return true if the site is not uniform around the block (it is handled individually or in a group)
         */
        abstract boolean isCandidate(int p);

        /**
         * Scores or draws the sites whose only events around the block are ambiguities.
         *
         * @return the log probability of their (drawn or existing) block states under the full conditional
         */
        abstract double sampleAmbiguitySites(Block block);

        void run(Block block) {
            otherChangedSites.clear();
            ambiguityLogProbability = sampleAmbiguitySites(block);
            uniformLogProbability = sampleUniformSites(block);
        }

        /**
         * Draws which sites of a group change, by thinning (every site in the group has the same probability
         * of changing), draws the new states of those, and scores the rest as unchanged.
         *
         * @param sites the group's sites by index, in pattern order (asked for in increasing order)
         * @return logProbability plus the log probability of the group's drawn states
         */
        double sampleGroup(Block block, UniformConditional uniform, int count, SiteLocator sites, double logProbability) {
            final double probabilityChanged = uniform.totalChangeMass;
            int changedCount = 0;
            if (probabilityChanged > 0.0) {
                final double logSkip = Math.log1p(-probabilityChanged);
                long index = -1;
                while (true) {
                    index = nextThinned(index, logSkip);
                    if (index >= count) {
                        break;
                    }
                    int p = sites.siteAt((int) index);
                    logProbability += blocks.sampleChangedUniformPattern(block, uniform, p);
                    otherChangedSites.add(p);
                    changedCount++;
                }
            }
            return logProbability + (count - changedCount) * uniform.logProbabilityUnchanged;
        }

        /**
         * Scores (or draws) the block states of the variable sites with no events around the block: each
         * keeps its state a throughout the block with probability P(unchanged | a). Draws use thinning, so
         * only the sites that change are visited (in expectation, as well as O(events)).
         *
         * @return the log probability of their (drawn or existing) block states under the full conditional
         */
        private double sampleUniformSites(Block block) {
            final double w = variableWeight;
            // by (rate class, state)
            final int classCount = uniformCounts.length;
            final UniformConditional[] uniform = new UniformConditional[classCount];
            final double[] probabilityChanged = new double[classCount];
            double maxProbabilityChanged = 0.0;
            for (int i = 0; i < classCount; i++) {
                if (uniformCounts[i] > 0) {
                    int rateClass = i / stateCount;
                    int a = i % stateCount;
                    uniform[i] = new UniformConditional(blocks.uniformStateConditional(block, a, w, blocks.rateOfClass(rateClass), rateClass), a);
                    probabilityChanged[i] = uniform[i].totalChangeMass;
                    maxProbabilityChanged = Math.max(maxProbabilityChanged, probabilityChanged[i]);
                }
            }

            changedUniformCount = 0;
            final int[] changedCounts = new int[classCount];
            double logProbability = 0.0;
            if (block.sample && maxProbabilityChanged > 0.0) {
                final double logSkip = Math.log1p(-maxProbabilityChanged);
                long i = -1;
                while (true) {
                    i = nextThinned(i, logSkip);
                    if (i >= variableSites.length) {
                        break;
                    }
                    int p = variableSites[(int) i];
                    if (isCandidate(p)) {
                        // a site with events: handled individually or in groups
                        continue;
                    }
                    int a = store.getState(referenceNode, p);
                    if (a < 0) {
                        continue;
                    }
                    int j = countIndex(p, a);
                    if (MathUtils.nextDouble() * maxProbabilityChanged >= probabilityChanged[j]) {
                        continue;
                    }
                    logProbability += blocks.sampleChangedUniformPattern(block, uniform[j], p);
                    changedCounts[j]++;
                    changedUniformCount++;
                    otherChangedSites.add(p);
                }
            }
            for (int i = 0; i < classCount; i++) {
                if (uniformCounts[i] > 0) {
                    logProbability += (uniformCounts[i] - changedCounts[i]) * uniform[i].logProbabilityUnchanged;
                }
            }
            return logProbability;
        }
    }

    /**
     * The sites of a group, by index in pattern order.
     */
    private interface SiteLocator {
        int siteAt(int index);
    }

    /**
     * @return the index of the next site visited by thinning, when each is skipped with probability
     * exp(logSkip): the number of sites skipped before it is geometric
     */
    private static long nextThinned(long index, double logSkip) {
        double u = 1.0 - MathUtils.nextDouble(); // in (0, 1]
        return index + 1 + (logSkip == Double.NEGATIVE_INFINITY ? 0L : (long) Math.floor(Math.log(u) / logSkip));
    }

    /**
     * The key of an ambiguity group: which of the block's neighbours are unambiguous, and the state.
     */
    private static final class AmbiguityKey {
        long[] present;
        int state;
        int hash;

        void set(long[] present, int state) {
            this.present = present;
            this.state = state;
            this.hash = 31 * Arrays.hashCode(present) + state;
        }

        @Override
        public int hashCode() {
            return hash;
        }

        @Override
        public boolean equals(Object o) {
            AmbiguityKey other = (AmbiguityKey) o;
            return state == other.state && Arrays.equals(present, other.present);
        }
    }

    /**
     * The sparse pass of a path move: its candidate sites (those with events around the block) are stamped in
     * siteStamps, and the ambiguity-only ones are listed, in pattern order.
     */
    private final class PathPass extends SparsePass {
        final int[] substitutionSites; // the candidates handled individually
        final int[] ambiguitySites;
        final int candidateStamp;

        PathPass(int referenceNode, int[] uniformCounts, int[] substitutionSites, int[] ambiguitySites,
                 int candidateStamp) {
            super(referenceNode, uniformCounts);
            this.substitutionSites = substitutionSites;
            this.ambiguitySites = ambiguitySites;
            this.candidateStamp = candidateStamp;
        }

        @Override
        boolean isCandidate(int p) {
            return siteStamps[p] == candidateStamp;
        }

        /**
         * Every unambiguous node around and in the block has the same state a at an ambiguity-only site, and
         * its block states stay at a (where the block has data below it) with a probability that depends only
         * on a and on which of the block's neighbours are ambiguous, so the sites are handled in groups as the
         * uniform sites are.
         */
        @Override
        double sampleAmbiguitySites(Block block) {
            if (ambiguitySites.length == 0) {
                return 0.0;
            }
            final int blockSize = block.nodes.length;
            final byte[][] blockStates = block.states;
            final int[] neighbours = block.neighbours;
            final int neighbourCount = neighbours.length;
            final byte[][] neighbourStates = new byte[neighbourCount][];
            for (int j = 0; j < neighbourCount; j++) {
                neighbourStates[j] = store.getStates(neighbours[j]);
            }

            // Group the sites. Which neighbours are ambiguous only changes at the ends of their ambiguous runs,
            // so the sites (in pattern order) are swept against those, and a site's state is read from one
            // unambiguous neighbour (every unambiguous node around the block has the same state).
            final int words = (neighbourCount + 63) / 64;
            final long[] events = getRunEnds(neighbours);
            final int eventCount = events.length;

            final Map<AmbiguityKey, AmbiguityGroup> groups = new HashMap<AmbiguityKey, AmbiguityGroup>();
            final AmbiguityKey probe = new AmbiguityKey();
            final long[] present = new long[words];
            for (int j = 0; j < neighbourCount; j++) {
                present[j >> 6] |= 1L << (j & 63);
            }
            int nextEvent = 0;
            boolean presenceChanged = true;
            int firstPresent = -1;
            // the groups of the current presence, by (rate class, state)
            final AmbiguityGroup[] groupsByState = new AmbiguityGroup[rateClassCount * stateCount];
            double logProbability = 0.0;
            for (int p : ambiguitySites) {
                while (nextEvent < eventCount && (int) (events[nextEvent] >>> 32) <= p) {
                    // the start of an ambiguous run makes the neighbour ambiguous, the end unambiguous again
                    int j = (int) events[nextEvent];
                    present[j >> 6] ^= 1L << (j & 63);
                    nextEvent++;
                    presenceChanged = true;
                }
                if (presenceChanged) {
                    firstPresent = -1;
                    for (int word = 0; word < words; word++) {
                        if (present[word] != 0L) {
                            firstPresent = (word << 6) + Long.numberOfTrailingZeros(present[word]);
                            break;
                        }
                    }
                    Arrays.fill(groupsByState, null);
                    presenceChanged = false;
                }
                int a;
                if (firstPresent >= 0) {
                    a = neighbourStates[firstPresent][p];
                } else {
                    // no data around the block (block nodes may still have states from before the move)
                    a = -1;
                    for (int i = 0; i < blockSize; i++) {
                        a = Math.max(a, blockStates[i][p]);
                    }
                }
                if (a < 0) {
                    // ambiguous throughout
                    setAll(block, p, -1, null);
                    continue;
                }
                final int i = countIndex(p, a);
                AmbiguityGroup current = groupsByState[i];
                if (current == null) {
                    probe.set(present, i);
                    current = groups.get(probe);
                    if (current == null) {
                        AmbiguityKey key = new AmbiguityKey();
                        key.set(present.clone(), i);
                        current = new AmbiguityGroup(a, p);
                        groups.put(key, current);
                    }
                    groupsByState[i] = current;
                }
                current.sites.add(p);
            }

            final double w = variableWeight;
            for (AmbiguityGroup group : groups.values()) {
                final int count = group.sites.count;
                final Conditional conditional = blocks.siteConditional(block, group.representative, w,
                        delegate.getPatternRate(group.representative));
                if (!conditional.present[block.topIndex]) {
                    // no data below the block
                    for (int i = 0; i < count; i++) {
                        setAll(block, group.sites.values[i], -1, null);
                    }
                    continue;
                }
                final UniformConditional uniform = new UniformConditional(conditional, group.state);
                if (!block.sample) {
                    logProbability += count * uniform.logProbabilityUnchanged;
                    continue;
                }
                // the unchanged configuration: unambiguous nodes take a (which they may not have had before the move)
                for (int i = 0; i < count; i++) {
                    setAll(block, group.sites.values[i], group.state, conditional.present);
                }
                logProbability = sampleGroup(block, uniform, count, group, logProbability);
            }
            return logProbability;
        }

        /**
         * Sets a site's block states to the state (where present, or everywhere if present is null), noting
         * the site as changed if any of them change. Does nothing if not sampling.
         */
        private void setAll(Block block, int p, int state, boolean[] present) {
            if (!block.sample) {
                return;
            }
            final byte[][] states = block.states;
            boolean changed = false;
            for (int i = 0; i < states.length; i++) {
                byte value = (byte) (present == null || present[i] ? state : -1);
                if (states[i][p] != value) {
                    block.setState(i, p, value);
                    changed = true;
                }
            }
            if (changed) {
                otherChangedSites.add(p);
            }
        }
    }

    /**
     * The ambiguity-only sites of a path move that share a state and a pattern of ambiguous neighbours (and
     * so a conditional).
     */
    private static final class AmbiguityGroup implements SiteLocator {
        final int state;
        final int representative;
        final IntBuffer sites = new IntBuffer();

        AmbiguityGroup(int state, int representative) {
            this.state = state;
            this.representative = representative;
        }

        @Override
        public int siteAt(int index) {
            return sites.values[index];
        }
    }

    /**
     * @return the ends of the nodes' runs of ambiguous patterns, in pattern order, as (position << 32) | j
     * for the j-th node
     */
    private long[] getRunEnds(int[] nodes) {
        int eventCount = 0;
        for (int n : nodes) {
            eventCount += getAmbiguousRuns(n).length;
        }
        final long[] events = new long[eventCount];
        int k = 0;
        for (int j = 0; j < nodes.length; j++) {
            for (int position : getAmbiguousRuns(nodes[j])) {
                events[k++] = ((long) position << 32) | j;
            }
        }
        Arrays.sort(events);
        return events;
    }

    /**
     * The sparse pass of a local move, for the block {P}: the sites at which some of P's neighbours are
     * ambiguous (and that are not substitution or hot sites) are swept as segments over which the same
     * neighbours are ambiguous, and grouped by their state, without lists of sites. P's states only need
     * writing where its presence changes, and the sites that change otherwise are located when drawn. The
     * groups, and the order of their sites, are those a path pass would form.
     */
    private final class LocalPass extends SparsePass {
        // the sites not handled here: constant, substitution (handled individually) and hot (handled as chains)
        final long[] excluded;
        byte[][] neighbourStates;      // P's neighbours' states (set when the ambiguity-only sites are swept)

        LocalPass(int referenceNode, int[] uniformCounts, int[] substitutionSites, int[] hotSites) {
            super(referenceNode, uniformCounts);
            this.excluded = constantBits.clone();
            for (int p : substitutionSites) {
                excluded[p >>> 6] |= 1L << p;
            }
            for (int p : hotSites) {
                excluded[p >>> 6] |= 1L << p;
            }
        }

        /**
         * @return true for a substitution site, a hot site, or one at which one of P's neighbours is ambiguous
         */
        @Override
        boolean isCandidate(int p) {
            if (isExcluded(p)) {
                return true;
            }
            for (byte[] states : neighbourStates) {
                if (states[p] < 0) {
                    return true;
                }
            }
            return false;
        }

        private boolean isExcluded(int p) {
            return ((excluded[p >>> 6] >>> p) & 1L) != 0L;
        }

        /**
         * Adds the number of sites of [start, end) in each (rate class, state) (all must be unambiguous) that are
         * not excluded.
         */
        private void countStates(byte[] states, int start, int end, int[] counts) {
            if (rateClassCount > 1) {
                final byte[] categories = store.getCategories();
                for (int p = start; p < end; p++) {
                    if (!isExcluded(p)) {
                        counts[categories[p] * stateCount + states[p]]++;
                    }
                }
                return;
            }
            int p = start;
            while (p < end) {
                final long bits = excluded[p >>> 6];
                final int limit = Math.min(end, ((p >>> 6) + 1) << 6);
                if (bits == 0L) {
                    for (; p < limit; p++) {
                        counts[states[p]]++;
                    }
                } else {
                    for (; p < limit; p++) {
                        if (((bits >>> p) & 1L) == 0L) {
                            counts[states[p]]++;
                        }
                    }
                }
            }
        }

        @Override
        double sampleAmbiguitySites(Block block) {
            // P's neighbours: its parent (if any), then its children
            final int[] neighbours = block.neighbours;
            final int neighbourCount = neighbours.length;
            int childBits = 0;
            int xBit = 0;
            for (int k = block.topParent == null ? 0 : 1; k < neighbourCount; k++) {
                childBits |= 1 << k;
                if (neighbours[k] == referenceNode) {
                    xBit = 1 << k;
                }
            }
            neighbourStates = new byte[neighbourCount][];
            for (int j = 0; j < neighbourCount; j++) {
                neighbourStates[j] = store.getStates(neighbours[j]);
            }

            // the ends of the neighbours' ambiguous runs, in pattern order
            final long[] events = getRunEnds(neighbours);
            final int eventCount = events.length;
            // where P has data before this pass's writes
            final int[] pRuns = getAmbiguousRuns(block.nodes[0]);

            final int fullMask = (1 << neighbourCount) - 1;
            final Map<AmbiguityKey, LocalAmbiguityGroup> groups = new HashMap<AmbiguityKey, LocalAmbiguityGroup>();
            final AmbiguityKey probe = new AmbiguityKey();
            final long[] present = new long[1];
            // by (rate class, state)
            final int classCount = rateClassCount * stateCount;
            final int[] segmentCounts = new int[classCount];
            final int[] firstSites = new int[classCount];
            final int[] newStates = new int[classCount];
            int mask = fullMask;
            int start = 0;
            int nextEvent = 0;
            int pRun = 0;
            while (start < patternCount) {
                int end = nextEvent < eventCount ? (int) (events[nextEvent] >>> 32) : patternCount;
                if (end > start && mask != fullMask) {
                    // a segment over which the same neighbours are ambiguous
                    int first = Integer.numberOfTrailingZeros(mask);
                    boolean pPresent = (mask & childBits) != 0;
                    if (first < neighbourCount) {
                        present[0] = mask;
                        final byte[] states = neighbourStates[first];
                        Arrays.fill(segmentCounts, 0);
                        countStates(states, start, end, segmentCounts);
                        // the groups of the segment's states, in the order of their first sites (the order in
                        // which new groups are created, as when the sites were taken one at a time)
                        int stateCountInSegment = 0;
                        Arrays.fill(firstSites, -1);
                        for (int i = 0; i < classCount; i++) {
                            if (segmentCounts[i] > 0) {
                                stateCountInSegment++;
                            }
                        }
                        int found = 0;
                        for (int p = start; found < stateCountInSegment; p++) {
                            if (!isExcluded(p)) {
                                int i = countIndex(p, states[p]);
                                if (firstSites[i] < 0) {
                                    firstSites[i] = p;
                                    newStates[found++] = i;
                                }
                            }
                        }
                        for (int j = 0; j < found; j++) {
                            final int i = newStates[j];
                            probe.set(present, i);
                            LocalAmbiguityGroup group = groups.get(probe);
                            if (group == null) {
                                AmbiguityKey key = new AmbiguityKey();
                                key.set(present.clone(), i);
                                group = new LocalAmbiguityGroup(i % stateCount, i / stateCount, firstSites[i], first);
                                groups.put(key, group);
                            }
                            group.addSegment(start, end);
                            group.count += segmentCounts[i];
                        }
                    }
                    if (block.sample) {
                        // P's states where its presence changes (elsewhere it already has the neighbours' state)
                        while (pRun < pRuns.length && pRuns[pRun + 1] <= start) {
                            pRun += 2;
                        }
                        int position = start;
                        for (int r = pRun; position < end; r += 2) {
                            int runStart = r < pRuns.length ? Math.max(pRuns[r], start) : end;
                            int runEnd = r < pRuns.length ? Math.min(pRuns[r + 1], end) : end;
                            // [position, runStart): P had data; [runStart, runEnd): it had none
                            if (!pPresent) {
                                for (int p = position; p < Math.min(runStart, end); p++) {
                                    if (!isExcluded(p)) {
                                        block.setState(0, p, -1);
                                        otherChangedSites.add(p);
                                    }
                                }
                            } else {
                                for (int p = Math.max(runStart, position); p < runEnd; p++) {
                                    if (!isExcluded(p)) {
                                        block.setState(0, p, neighbourStates[first][p]);
                                        otherChangedSites.add(p);
                                    }
                                }
                            }
                            position = Math.max(runEnd, position);
                            if (r >= pRuns.length) {
                                break;
                            }
                        }
                    }
                }
                if (nextEvent >= eventCount) {
                    break;
                }
                // the start of a run makes a neighbour ambiguous, its end unambiguous again
                mask ^= 1 << (int) events[nextEvent];
                nextEvent++;
                start = end;
            }

            final double w = variableWeight;
            double logProbability = 0.0;
            for (Map.Entry<AmbiguityKey, LocalAmbiguityGroup> entry : groups.entrySet()) {
                final LocalAmbiguityGroup group = entry.getValue();
                final int groupMask = (int) entry.getKey().present[0];
                // these sites are not uniform: take them out of the counts of the uniform sites (read at X)
                if ((groupMask & xBit) != 0) {
                    uniformCounts[group.rateClass * stateCount + group.state] -= group.count;
                }
                final Conditional conditional = blocks.siteConditional(block, group.representative, w,
                        delegate.getPatternRate(group.representative));
                if (!conditional.present[0]) {
                    // no data below P: its states were cleared above
                    continue;
                }
                final UniformConditional uniform = new UniformConditional(conditional, group.state);
                if (!block.sample) {
                    logProbability += group.count * uniform.logProbabilityUnchanged;
                    continue;
                }
                group.resetCursor();
                logProbability = sampleGroup(block, uniform, group.count, index -> locate(group, index), logProbability);
            }
            return logProbability;
        }

        /**
         * @return the site of the group with the given index (in pattern order); indices must be asked for in
         * increasing order after resetCursor
         */
        private int locate(LocalAmbiguityGroup group, int index) {
            final byte[] states = neighbourStates[group.firstPresent];
            final byte[] categories = store.getCategories();
            while (true) {
                int segmentStart = group.segments.values[2 * group.cursorSegment];
                int segmentEnd = group.segments.values[2 * group.cursorSegment + 1];
                if (group.cursorSite < segmentStart) {
                    group.cursorSite = segmentStart;
                }
                for (; group.cursorSite < segmentEnd; group.cursorSite++) {
                    int p = group.cursorSite;
                    if (states[p] == group.state && (categories == null || categories[p] == group.rateClass) &&
                            !isExcluded(p)) {
                        if (group.cursorIndex == index) {
                            group.cursorIndex++;
                            group.cursorSite++;
                            return p;
                        }
                        group.cursorIndex++;
                    }
                }
                group.cursorSegment++;
            }
        }
    }

    /**
     * The sites of a local move's segments that share a pattern of ambiguous neighbours and a state.
     */
    private static final class LocalAmbiguityGroup {
        final int state;
        final int rateClass;
        final int representative;
        final int firstPresent;   // the neighbour the state is read from
        final IntBuffer segments = new IntBuffer(); // [start, end) of each segment with sites in the group
        int count = 0;
        int cursorSegment;
        int cursorSite;
        int cursorIndex;

        LocalAmbiguityGroup(int state, int rateClass, int representative, int firstPresent) {
            this.state = state;
            this.rateClass = rateClass;
            this.representative = representative;
            this.firstPresent = firstPresent;
        }

        void addSegment(int start, int end) {
            segments.add(start);
            segments.add(end);
        }

        void resetCursor() {
            cursorSegment = 0;
            cursorSite = 0;
            cursorIndex = 0;
        }
    }
}
