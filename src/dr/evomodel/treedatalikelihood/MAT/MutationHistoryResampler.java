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
 * block states stay unchanged with a probability that depends only on a. So the sites that change are
 * drawn by thinning (stepping through the variable sites with geometric jumps at the highest probability
 * of a change, and keeping each site visited with its own probability over that), and the unchanged ones
 * contribute n_a log P(unchanged | a) for the number n_a of them in each state. Sites whose only events
 * around the block are ambiguities are handled in the same way, grouped by their state and by which of
 * the block's neighbours are ambiguous there: there are typically many such sites (runs of N) but few
 * groups. A move then costs O(events around the block + changes), not O(sites).
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
    private final BranchHistoryModel historyModel;
    private final Tree tree;
    private final int stateCount;
    private final int patternCount;
    private final boolean localBlock;

    // patterns grouped by weight, so conditionals can be looked up by [state][weight class]
    private final int[] weightClass;
    private final double[] classWeights;
    private final int[] variableSites;
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
    private final double[] newSiteLogProbabilities;
    private final boolean[] resampledSites;
    private final int[] changedSiteStamps;
    private int changedSiteStamp = 0;
    private long[] siteBits;
    private int[] oldMarks;
    private int[] newMarks;
    private int[] hotMarks;
    private int localStamp = 0;

    // per-node work arrays and caches
    private final boolean[] inBlock;
    private int[] nodeBlockIndex;
    private int[] nodeMarks;           // a node is marked if its entry equals nodeStamp
    private int nodeStamp = 0;
    private int[][] ambiguousRuns;
    private long[] ambiguousRunsVersion;
    private int[][] stateCountsCache;
    private long[] stateCountsVersion;

    private int chosenState;           // left by choose

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
        this.historyModel = delegate.getHistoryModel();
        this.tree = delegate.getTree();
        this.stateCount = delegate.getPatternList().getDataType().getStateCount();
        this.patternCount = delegate.getPatternList().getPatternCount();
        this.inBlock = new boolean[tree.getNodeCount()];
        this.siteClasses = new byte[patternCount];
        this.siteStamps = new int[patternCount];
        this.oldSiteLogProbabilities = new double[patternCount];
        this.newSiteLogProbabilities = new double[patternCount];
        this.resampledSites = new boolean[patternCount];
        this.changedSiteStamps = new int[patternCount];

        double[] weights = delegate.getPatternWeights();
        Map<Double, Integer> classes = new HashMap<Double, Integer>();
        this.weightClass = new int[patternCount];
        for (int p = 0; p < patternCount; p++) {
            Integer c = classes.get(weights[p]);
            if (c == null) {
                c = classes.size();
                classes.put(weights[p], c);
            }
            weightClass[p] = c;
        }
        this.classWeights = new double[classes.size()];
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
        this.localBlock = localBlock;
        if (localBlock && !sparse) {
            throw new IllegalArgumentException("The local block needs a branch history model without explicit " +
                    "hits, no per-site rates and equal pattern weights");
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

    @Override
    public void beforeSubtreeMove(NodeRef node, NodeRef sibling, NodeRef destination) {
        if (localBlock) {
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
            // the sites with events around the block; the others are handled together (see SparsePass)
            candidateSites = getCandidateSites(block);
            final DefaultMutationStore store = delegate.getMutableMutationStore();
            final int x = node.getNumber();
            final int[] uniformCounts = getStateCounts(x).clone();
            int substitutionCount = 0;
            int ambiguityCount = 0;
            int[] substitution = new int[candidateSites.length];
            int[] ambiguity = new int[candidateSites.length];
            for (int p : candidateSites) {
                if (delegate.isConstantPattern(p)) {
                    continue;
                }
                int a = store.getState(x, p);
                if (a >= 0) {
                    uniformCounts[a]--;
                }
                if (siteClasses[p] == AMBIGUITIES_ONLY) {
                    ambiguity[ambiguityCount++] = p;
                } else {
                    substitution[substitutionCount++] = p;
                }
            }
            // the ambiguity-only sites in pattern order, so that they can be swept against the neighbours' runs
            oldPathPass = new PathPass(x, uniformCounts, Arrays.copyOf(substitution, substitutionCount),
                    sortSites(ambiguity, ambiguityCount), stamp);
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

    private double afterPathMove() {
        if (block == null) {
            throw new IllegalStateException("afterSubtreeMove called without beforeSubtreeMove");
        }
        final int[] currentBlock = block;
        final PathPass oldPass = oldPathPass;
        block = null;
        oldPathPass = null;

        final DefaultMutationStore store = delegate.getMutableMutationStore();
        final List<Integer> reparented = new ArrayList<Integer>();
        IntBuffer synced = syncParents(movedNodes);
        for (int i = 0; i < synced.count; i++) {
            reparented.add(synced.values[i]);
        }

        // every candidate is resampled (the others are handled together when sampling sparsely)
        final int[] resampled = oldPass != null ? oldPass.substitutionSites : candidateSites;
        for (int p : resampled) {
            resampledSites[p] = true;
        }

        double logProposalOfOldStates = oldPass != null ?
                oldPass.uniformLogProbability + oldPass.ambiguityLogProbability : 0.0;
        for (int p : resampled) {
            logProposalOfOldStates += oldSiteLogProbabilities[p];
        }
        final PathPass newPass = oldPass == null ? null :
                new PathPass(oldPass.referenceNode, oldPass.uniformCounts, oldPass.substitutionSites,
                        oldPass.ambiguitySites, oldPass.candidateStamp);
        double logProposalOfNewStates = sampleBlock(currentBlock, true, siteClasses, resampled,
                newSiteLogProbabilities, newPass);
        lastChangedUniformCount = newPass == null ? 0 : newPass.changedUniformCount;
        delegate.setHistoryResampled();

        int[] changed = resampled;
        if (newPass != null && newPass.changedSites.count > 0) {
            // the other sites whose block states changed also need their branches updated
            changed = Arrays.copyOf(resampled, resampled.length + newPass.changedSites.count);
            int k = resampled.length;
            for (int i = 0; i < newPass.changedSites.count; i++) {
                changed[k++] = newPass.changedSites.values[i];
                resampledSites[newPass.changedSites.values[i]] = true;
            }
        }

        // the Fitch sets no longer match the topology; parsimony will recompute them if it is ever needed
        store.setFitchSetsStale(true);

        if (historyModel.hasExplicitHistories()) {
            // draw the hits of the resampled sites on every branch touching the block, given their end states
            Set<Integer> branches = new HashSet<Integer>(reparented);
            for (int n : getIncidentBranches(currentBlock)) {
                branches.add(n);
            }
            double logProposalOfNewHistories = sampleHistories(branches, resampled);
            for (int p : resampled) {
                resampledSites[p] = false;
            }
            candidateSites = null;
            delegate.markTopologySynchronised();
            return logProposalOfOldStates - logProposalOfNewStates + logProposalOfOldHistories - logProposalOfNewHistories;
        }

        // branches whose parent changed are recomputed in full; the others around the block changed only at the resampled sites
        changed = sortSites(changed, changed.length);
        Set<Integer> fullUpdates = new HashSet<Integer>(reparented);
        for (int n : fullUpdates) {
            store.updateBranch(n);
        }
        for (int n : currentBlock) {
            if (!fullUpdates.contains(n)) {
                store.updateBranchSites(n, resampledSites, changed);
            }
            NodeRef node = tree.getNode(n);
            for (int i = 0; i < tree.getChildCount(node); i++) {
                int c = tree.getChild(node, i).getNumber();
                if (!fullUpdates.contains(c) && !contains(currentBlock, c)) {
                    store.updateBranchSites(c, resampledSites, changed);
                }
            }
        }
        for (int p : changed) {
            resampledSites[p] = false;
        }
        candidateSites = null;
        delegate.markTopologySynchronised();

        return logProposalOfOldStates - logProposalOfNewStates;
    }

    /**
     * @return the regrafted parent plus the internal nodes on the path between a and b
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
        MutationStore store = delegate.getMutationStore();
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
        MutationStore store = delegate.getMutationStore();
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
        int x;
        int p;
        int s;
        int d;
        int newParent;                 // -1 if P becomes the root
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
        final DefaultMutationStore store = delegate.getMutableMutationStore();
        if (oldMarks == null) {
            oldMarks = new int[patternCount];
            newMarks = new int[patternCount];
            hotMarks = new int[patternCount];
        }
        localStamp++;
        final int markStamp = localStamp;

        final LocalMove move = new LocalMove();
        final NodeRef pNode = tree.getParent(node);
        final NodeRef gNode = tree.getParent(pNode);
        final NodeRef dParent = tree.getParent(destination);
        // P's new parent: D's parent, unless D is S (which then stays below G)
        final NodeRef newParentNode = dParent != null && dParent.getNumber() == pNode.getNumber() ? gNode : dParent;
        move.x = node.getNumber();
        move.p = pNode.getNumber();
        move.s = sibling.getNumber();
        move.d = destination.getNumber();
        move.newParent = newParentNode == null ? -1 : newParentNode.getNumber();
        move.presentX = SiteRuns.complement(getAmbiguousRuns(move.x), patternCount);

        // The old chain: G and each ancestor above it lose their states where X was their only data, up to
        // the first node that is above X's new position too.
        final int newAncestorStamp = nextNodeStamp();
        for (NodeRef n = newParentNode; n != null; n = tree.getParent(n)) {
            nodeMarks[n.getNumber()] = newAncestorStamp;
        }
        final IntBuffer chain = new IntBuffer();
        final List<int[]> levels = new ArrayList<int[]>();
        int[] hot = SiteRuns.intersect(move.presentX, getAmbiguousRuns(move.s));
        int previous = move.p;
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
        for (int site : move.oldHotSites) {
            hotMarks[site] = markStamp;
        }

        // The old pass: P with neighbours G, S and X.
        // (the sites at which a neighbour is ambiguous are found from their runs when P is scored or drawn)
        final IntBuffer oldSubstitutions = new IntBuffer();
        collectSubstitutionSites(store, gNode == null ? new int[]{move.x, move.s} : new int[]{move.x, move.s, move.p},
                oldMarks, markStamp, oldSubstitutions);
        final int[] oldSubstitutionSites = oldSubstitutions.toArray();
        final int[] oldUniformCounts = getUniformCounts(move.x, oldSubstitutionSites, move.oldHotSites, oldMarks, markStamp);

        // The new pass: P with neighbours N, D and X. P's state there is affected by the substitutions
        // between X and N (in the old topology, where they are on the path between them) and between D and
        // N, and by the runs of ambiguous patterns of X, D and N (after the old chain loses its states, since
        // D may be in it).
        final NodeRef anchor = newParentNode != null ? newParentNode : destination;
        final int xAncestorStamp = nextNodeStamp();
        for (NodeRef n = node; n != null; n = tree.getParent(n)) {
            nodeMarks[n.getNumber()] = xAncestorStamp;
        }
        final IntBuffer pathBranches = new IntBuffer();
        NodeRef mrca = anchor;
        while (nodeMarks[mrca.getNumber()] != xAncestorStamp) {
            pathBranches.add(mrca.getNumber());
            mrca = tree.getParent(mrca);
        }
        for (NodeRef n = node; n.getNumber() != mrca.getNumber(); n = tree.getParent(n)) {
            pathBranches.add(n.getNumber());
        }
        if (newParentNode != null) {
            pathBranches.add(move.d);
        }
        final IntBuffer newSubstitutions = new IntBuffer();
        collectSubstitutionSites(store, pathBranches.toArray(), newMarks, markStamp, newSubstitutions);
        move.newSubstitutionSites = newSubstitutions.toArray();
        move.newStamp = markStamp;

        // score P's current states (and the old chain's) under the old topology
        final LocalPass pass = new LocalPass(move.x, oldUniformCounts, oldMarks, markStamp, markStamp);
        final int[] substitutionSites = withoutMarked(oldSubstitutionSites, hotMarks, markStamp);
        double logProbability = sampleBlock(new int[]{move.p}, false, null, substitutionSites, null, pass);
        logProbability += sampleChains(move.p, move.oldChain, move.oldChainLevels, false);
        move.logProposalOfOldStates = logProbability;
        localMove = move;
    }

    private double afterLocalMove() {
        final LocalMove move = localMove;
        if (move == null) {
            throw new IllegalStateException("afterSubtreeMove called without beforeSubtreeMove");
        }
        localMove = null;
        final DefaultMutationStore store = delegate.getMutableMutationStore();
        final IntBuffer reparented = syncParents(new int[]{move.s, move.p, move.d});

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
        final LocalPass pass = new LocalPass(move.x,
                getUniformCounts(move.x, move.newSubstitutionSites, newHotSites, newMarks, move.newStamp),
                newMarks, move.newStamp, newHotStamp);
        final int[] substitutionSites = withoutMarked(move.newSubstitutionSites, hotMarks, newHotStamp);
        double logProposalOfNewStates = sampleBlock(new int[]{move.p}, true, null, substitutionSites, null, pass);
        logProposalOfNewStates += sampleChains(move.p, newChain, newChainLevels, true);
        lastChangedUniformCount = pass.changedUniformCount;
        delegate.setHistoryResampled();
        // the Fitch sets no longer match the topology; parsimony will recompute them if it is ever needed
        store.setFitchSetsStale(true);

        // the sites at which any state changed
        final IntBuffer changed = new IntBuffer();
        addFlagged(changed, substitutionSites, substitutionSites.length);
        addFlagged(changed, pass.changedSites.values, pass.changedSites.count);
        addFlagged(changed, move.oldHotSites, move.oldHotSites.length);
        addFlagged(changed, newHotSites, newHotSites.length);
        final int[] changedSites = sortSites(changed.values, changed.count);

        // the branches whose parent changed are recomputed in full; those whose ends changed only at the
        // changed sites are updated there
        final int touchedStamp = nextNodeStamp();
        for (int i = 0; i < reparented.count; i++) {
            store.updateBranch(reparented.values[i]);
            nodeMarks[reparented.values[i]] = touchedStamp;
        }
        // the other branches below the nodes whose states changed (and above them)
        final IntBuffer touched = new IntBuffer();
        addUnmarked(touched, move.x, touchedStamp);
        addWithChildren(touched, move.p, touchedStamp);
        for (int n : move.oldChain) {
            addWithChildren(touched, n, touchedStamp);
        }
        for (int n : newChain) {
            addWithChildren(touched, n, touchedStamp);
        }
        for (int i = 0; i < touched.count; i++) {
            store.updateBranchSites(touched.values[i], resampledSites, changedSites);
        }
        for (int site : changedSites) {
            resampledSites[site] = false;
        }
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
    private void collectSubstitutionSites(DefaultMutationStore store, int[] branches, int[] marks, int markStamp,
                                          IntBuffer out) {
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
        final DefaultMutationStore store = delegate.getMutableMutationStore();
        final int[] counts = getStateCounts(x).clone();
        for (int site : substitutions) {
            int a = store.getState(x, site);
            if (a >= 0) {
                counts[a]--;
            }
        }
        for (int site : hot) {
            int a = store.getState(x, site);
            if (a >= 0 && marks[site] != markStamp) {
                counts[a]--;
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

    private void addFlagged(IntBuffer buffer, int[] sites, int count) {
        for (int i = 0; i < count; i++) {
            if (!resampledSites[sites[i]]) {
                resampledSites[sites[i]] = true;
                buffer.add(sites[i]);
            }
        }
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
        if (nodeMarks == null) {
            nodeMarks = new int[tree.getNodeCount()];
        }
        return ++nodeStamp;
    }

    /**
     * Gives the store the tree's parents of the nodes a subtree move reparents (the only ones it changes).
     *
     * @return the nodes whose parent changed
     */
    private IntBuffer syncParents(int[] nodes) {
        final DefaultMutationStore store = delegate.getMutableMutationStore();
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
        final DefaultMutationStore store = delegate.getMutableMutationStore();
        int[] internal = getInternalNodes();
        sampleBlock(internal, true);
        if (historyModel.hasExplicitHistories()) {
            int[] sites = getVariableSites(internal);
            for (int p : sites) {
                resampledSites[p] = true;
            }
            List<Integer> branches = new ArrayList<Integer>();
            for (int n : getNonRootNodes()) {
                branches.add(n);
            }
            sampleHistories(branches, sites);
            for (int p : sites) {
                resampledSites[p] = false;
            }
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
        final DefaultMutationStore store = delegate.getMutableMutationStore();
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
     * Draws the hits of the sites flagged in resampledSites (listed in sites) on the branches above the
     * given nodes, given their end states, keeping the events of the other sites.
     *
     * @return the log probability of the drawn hits
     */
    private double sampleHistories(Iterable<Integer> branches, int[] sites) {
        final DefaultMutationStore store = delegate.getMutableMutationStore();
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
                if (!resampledSites[old.getPattern(k)]) {
                    for (int i = k; i < siteEnd; i++) {
                        builder.addSubstitution(old.getPattern(i), old.getFromState(i), old.getToState(i),
                                old.hasTimes() ? old.getTime(i) : (i - k + 1.0) / (siteEnd - k + 1.0));
                    }
                }
                k = siteEnd;
            }
            for (int i = old.getSubstitutionCount(); i < old.size(); i++) {
                if (!resampledSites[old.getPattern(i)]) {
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
        final DefaultMutationStore store = delegate.getMutableMutationStore();
        if (ambiguousRuns == null) {
            ambiguousRuns = new int[tree.getNodeCount()][];
            ambiguousRunsVersion = new long[tree.getNodeCount()];
        }
        long version = store.getStateVersion(n);
        if (ambiguousRuns[n] == null || ambiguousRunsVersion[n] != version) {
            ambiguousRuns[n] = SiteRuns.ambiguous(store.getStates(n), patternCount);
            ambiguousRunsVersion[n] = version;
        }
        return ambiguousRuns[n];
    }

    /**
     * @return the number of the variable sites in each state at the node (cached for its current states;
     * callers must not modify it)
     */
    private int[] getStateCounts(int n) {
        final DefaultMutationStore store = delegate.getMutableMutationStore();
        if (stateCountsCache == null) {
            stateCountsCache = new int[tree.getNodeCount()][];
            stateCountsVersion = new long[tree.getNodeCount()];
            Arrays.fill(stateCountsVersion, -1L);
        }
        long version = store.getStateVersion(n);
        if (stateCountsCache[n] == null || stateCountsVersion[n] != version) {
            int[] counts = stateCountsCache[n] == null ? new int[stateCount] : stateCountsCache[n];
            Arrays.fill(counts, 0);
            byte[] states = store.getStates(n);
            for (int p : variableSites) {
                if (states[p] >= 0) {
                    counts[states[p]]++;
                }
            }
            stateCountsCache[n] = counts;
            stateCountsVersion[n] = version;
        }
        return stateCountsCache[n];
    }

    /**
     * @return the given sites in pattern order (each site at most once)
     */
    private int[] sortSites(int[] sites, int count) {
        if (siteBits == null) {
            siteBits = new long[(patternCount + 63) / 64];
        }
        for (int i = 0; i < count; i++) {
            siteBits[sites[i] >> 6] |= 1L << (sites[i] & 63);
        }
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
    // Blocks
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
     * As {@link #sampleBlock(int[], boolean, byte[])} but only for the given sites (all sites if null),
     * optionally recording each site's log probability, and then, if a sparse pass is given, for the
     * sites it handles (those with no events, or only ambiguities, around the block).
     */
    private double sampleBlock(int[] blockNodes, boolean sample, byte[] patternClasses,
                               int[] sites, double[] siteLogProbabilities, SparsePass pass) {
        delegate.ensureRates();
        final Block block = new Block(blockNodes, sample);
        final DefaultMutationStore store = block.store;
        final int blockSize = block.nodes.length;
        final int topIndex = block.topIndex;
        final NodeRef topParent = block.topParent;
        final int[][] children = block.children;
        final int[] parentIndex = block.parentIndex;
        final byte[][] blockStates = block.states;
        final double siteRateScale = delegate.getSiteRateScale();
        final double[] patternWeights = delegate.getPatternWeights();

        // The full conditional of a pattern's block states depends only on its weight, its rate and the
        // states adjacent to the block, so it is computed once for each distinct combination of these.
        // Draws are still made independently for each pattern.
        final Map<ConditionalKey, Conditional> cache = new HashMap<ConditionalKey, Conditional>();
        final Map<ConditionalKey, UniformConditional> uniformCache = new HashMap<ConditionalKey, UniformConditional>();
        // without per-site rates, conditionals for patterns with no events depend only on (state, weight)
        final UniformConditional[][] uniformByClass = delegate.hasSiteRates() ? null :
                new UniformConditional[stateCount][classWeights.length];
        int boundaryCount = 1;
        for (int i = 0; i < blockSize; i++) {
            for (int c : children[i]) {
                if (!inBlock[c]) {
                    boundaryCount++;
                }
            }
        }
        final byte[] boundary = new byte[boundaryCount];

        double logProbability = 0.0;

        final int siteCount = sites == null ? patternCount : sites.length;
        int lastSite = -1;
        double lastLogProbability = 0.0;
        for (int siteIndex = 0; siteIndex < siteCount; siteIndex++) {
            final int p = sites == null ? siteIndex : sites[siteIndex];
            if (siteLogProbabilities != null) {
                if (lastSite >= 0) {
                    siteLogProbabilities[lastSite] = logProbability - lastLogProbability;
                }
                lastSite = p;
                lastLogProbability = logProbability;
            }
            final double w = patternWeights[p];
            final double nu = delegate.getSiteRate(p, siteRateScale);

            if (delegate.isConstantPattern(p)) {
                // constant sites are assumed to have no mutations: their states are never changed
                continue;
            }
            final byte patternClass = patternClasses == null ? HAS_SUBSTITUTIONS : patternClasses[p];
            if (patternClass == NO_EVENTS) {
                // every node in and around the block has this state (states are not yet overwritten)
                int a = blockStates[0][p];
                if (a < 0) {
                    continue;
                }
                UniformConditional uniform;
                if (uniformByClass != null) {
                    uniform = uniformByClass[a][weightClass[p]];
                    if (uniform == null) {
                        uniform = new UniformConditional(new Conditional(block, p, w, nu), a);
                        uniformByClass[a][weightClass[p]] = uniform;
                    }
                } else {
                    ConditionalKey key = new ConditionalKey(new byte[]{(byte) a}, w, nu);
                    uniform = uniformCache.get(key);
                    if (uniform == null) {
                        uniform = new UniformConditional(new Conditional(block, p, w, nu), a);
                        uniformCache.put(key, uniform);
                    }
                }
                if (!sample || MathUtils.nextDouble() < uniform.probabilityUnchanged) {
                    logProbability += uniform.logProbabilityUnchanged;
                } else {
                    logProbability += sampleChangedUniformPattern(block, uniform, p);
                }
                continue;
            }

            int k = 0;
            boundary[k++] = topParent == null ? -2 : (byte) store.getState(topParent.getNumber(), p);
            for (int i = 0; i < blockSize; i++) {
                for (int c : children[i]) {
                    if (!inBlock[c]) {
                        boundary[k++] = (byte) store.getState(c, p);
                    }
                }
            }
            ConditionalKey key = new ConditionalKey(boundary.clone(), w, nu);
            Conditional conditional = cache.get(key);
            if (conditional == null) {
                conditional = new Conditional(block, p, w, nu);
                cache.put(key, conditional);
            }

            if (patternClass == AMBIGUITIES_ONLY) {
                // every unambiguous node in and around the block has the same state, before and after the move
                int a = -1;
                for (byte b : boundary) {
                    a = Math.max(a, b);
                }
                for (int i = 0; i < blockSize; i++) {
                    a = Math.max(a, blockStates[i][p]);
                }
                if (a < 0 || !conditional.present[topIndex]) {
                    if (sample) {
                        for (int i = 0; i < blockSize; i++) {
                            block.setState(i, p, -1);
                        }
                    }
                    continue;
                }
                UniformConditional uniform = uniformCache.get(key);
                if (uniform == null) {
                    uniform = new UniformConditional(conditional, a);
                    uniformCache.put(key, uniform);
                }
                if (!sample) {
                    logProbability += uniform.logProbabilityUnchanged;
                    continue;
                }
                // the unchanged configuration: unambiguous nodes take a (which they may not have had before the move)
                for (int i = 0; i < blockSize; i++) {
                    block.setState(i, p, (byte) (conditional.present[i] ? a : -1));
                }
                if (MathUtils.nextDouble() < uniform.probabilityUnchanged) {
                    logProbability += uniform.logProbabilityUnchanged;
                } else {
                    logProbability += sampleChangedUniformPattern(block, uniform, p);
                }
                continue;
            }

            if (!conditional.present[topIndex]) {
                if (sample) {
                    for (int i = 0; i < blockSize; i++) {
                        block.setState(i, p, -1);
                    }
                }
                continue;
            }
            logProbability += choose(conditional.topLogP(), conditional.topP(), sample, blockStates[topIndex][p]);
            if (sample) {
                block.setState(topIndex, p, (byte) chosenState);
            }
            // pre-order: each child given its (already chosen) parent
            for (int i = topIndex - 1; i >= 0; i--) {
                if (!conditional.present[i]) {
                    if (sample) {
                        block.setState(i, p, -1);
                    }
                    continue;
                }
                int s = blockStates[parentIndex[i]][p];
                logProbability += choose(conditional.logRow(i, s), conditional.row(i, s), sample, blockStates[i][p]);
                if (sample) {
                    block.setState(i, p, (byte) chosenState);
                }
            }
        }

        if (siteLogProbabilities != null && lastSite >= 0) {
            siteLogProbabilities[lastSite] = logProbability - lastLogProbability;
        }

        if (pass != null) {
            pass.run(block, siteRateScale);
            logProbability += pass.ambiguityLogProbability;
            logProbability += pass.uniformLogProbability;
        }

        block.release();
        return logProbability;
    }

    /**
     * One call of sampleBlock: the block's nodes and branches in the current tree, the log branch densities
     * (for the rate of the patterns being handled), and the block's states, which are only changed through
     * {@link #setState} so that the store's node summaries stay up to date. The block's nodes are marked in
     * inBlock until {@link #release}.
     */
    private final class Block {
        final boolean sample;
        final DefaultMutationStore store;
        final int[] nodes;          // in post-order, the top last
        final int topIndex;
        final NodeRef topParent;
        final int[][] children;     // the children of each block node
        final int[][] childBranch;  // the index of the branch above each child
        final int[] index;          // node number -> block index (only read for block nodes)
        final int[] parentIndex;
        final int topBranch;        // the index of the branch above the top node
        final double[] branchLengths;
        final double[][] logF;      // the log branch densities, for the rate logFRate
        double logFRate = Double.NaN;
        final double[][][] linearMatrices;
        final double[] logFrequencies;
        final byte[][] states;

        Block(int[] blockNodes, boolean sample) {
            this.sample = sample;
            this.store = delegate.getMutableMutationStore();
            this.nodes = getBlockPostOrder(blockNodes);
            final int blockSize = nodes.length;
            this.topIndex = blockSize - 1;
            this.topParent = tree.getParent(tree.getNode(nodes[topIndex]));
            for (int n : nodes) {
                inBlock[n] = true;
            }

            this.children = new int[blockSize][];
            if (nodeBlockIndex == null) {
                nodeBlockIndex = new int[tree.getNodeCount()];
            }
            this.index = nodeBlockIndex;
            for (int i = 0; i < blockSize; i++) {
                index[nodes[i]] = i;
            }
            int branchTotal = 1;
            for (int i = 0; i < blockSize; i++) {
                branchTotal += tree.getChildCount(tree.getNode(nodes[i]));
            }
            this.branchLengths = new double[branchTotal];
            int branchCounter = 0;
            this.childBranch = new int[blockSize][];
            for (int i = 0; i < blockSize; i++) {
                NodeRef node = tree.getNode(nodes[i]);
                children[i] = new int[tree.getChildCount(node)];
                childBranch[i] = new int[children[i].length];
                for (int j = 0; j < children[i].length; j++) {
                    NodeRef child = tree.getChild(node, j);
                    children[i][j] = child.getNumber();
                    childBranch[i][j] = branchCounter;
                    branchLengths[branchCounter++] = delegate.getBranchLength(child);
                }
            }
            this.parentIndex = new int[blockSize];
            parentIndex[topIndex] = -1;
            for (int i = 0; i < blockSize; i++) {
                for (int c : children[i]) {
                    if (inBlock[c]) {
                        parentIndex[index[c]] = i;
                    }
                }
            }
            this.topBranch = branchCounter;
            branchLengths[branchCounter] = delegate.getBranchLength(tree.getNode(nodes[topIndex]));

            double[] frequencies = delegate.getRootFrequencyModel().getFrequencies();
            this.logFrequencies = new double[stateCount];
            for (int a = 0; a < stateCount; a++) {
                logFrequencies[a] = Math.log(frequencies[a]);
            }
            this.logF = new double[branchTotal][stateCount * stateCount];
            this.linearMatrices = new double[classWeights.length][][];
            ensureLogF(delegate.getSiteRateScale());

            this.states = new byte[blockSize][];
            for (int i = 0; i < blockSize; i++) {
                states[i] = sample ? store.getStatesForUpdate(nodes[i]) : store.getStates(nodes[i]);
            }
        }

        /**
         * Makes logF hold the log branch densities for the given relative rate.
         */
        void ensureLogF(double nu) {
            if (nu != logFRate) {
                for (int e = 0; e < logF.length; e++) {
                    historyModel.getLogTransitionDensities(branchLengths[e], nu, logF[e]);
                }
                logFRate = nu;
            }
        }

        /**
         * @return exp(w logF) for every branch, for the weight class's w (at the common rate, i.e. without
         * per-site rates)
         */
        double[][] getLinearMatrices(int weightClassIndex) {
            double[][] matrices = linearMatrices[weightClassIndex];
            if (matrices == null) {
                double w = classWeights[weightClassIndex];
                matrices = new double[logF.length][stateCount * stateCount];
                for (int e = 0; e < logF.length; e++) {
                    for (int j = 0; j < stateCount * stateCount; j++) {
                        matrices[e][j] = Math.exp(w * logF[e][j]);
                    }
                }
                linearMatrices[weightClassIndex] = matrices;
            }
            return matrices;
        }

        /**
         * Sets the state of block node i at a pattern, through the store.
         */
        void setState(int i, int p, int state) {
            if (states[i][p] != state) {
                store.setState(nodes[i], p, state);
            }
        }

        void release() {
            for (int n : nodes) {
                inBlock[n] = false;
            }
        }
    }

    /**
     * @return the block's nodes in post-order of the current tree, starting from its unique top node
     */
    private int[] getBlockPostOrder(int[] blockNodes) {
        for (int n : blockNodes) {
            inBlock[n] = true;
        }
        int top = -1;
        for (int n : blockNodes) {
            NodeRef parent = tree.getParent(tree.getNode(n));
            if (parent == null || !inBlock[parent.getNumber()]) {
                if (top != -1) {
                    clearBlock(blockNodes);
                    throw new IllegalStateException("Resampling block is not connected");
                }
                top = n;
            }
        }
        // pre-order by DFS restricted to the block, then reverse
        int[] order = new int[blockNodes.length];
        int[] stack = new int[blockNodes.length];
        int sp = 0;
        int k = blockNodes.length;
        stack[sp++] = top;
        while (sp > 0) {
            int n = stack[--sp];
            order[--k] = n;
            NodeRef node = tree.getNode(n);
            for (int i = 0; i < tree.getChildCount(node); i++) {
                int c = tree.getChild(node, i).getNumber();
                if (inBlock[c]) {
                    stack[sp++] = c;
                }
            }
        }
        clearBlock(blockNodes);
        if (k != 0) {
            throw new IllegalStateException("Resampling block is not connected");
        }
        return order;
    }

    private void clearBlock(int[] blockNodes) {
        for (int n : blockNodes) {
            inBlock[n] = false;
        }
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
        final IntBuffer changedSites = new IntBuffer();
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
        abstract double sampleAmbiguitySites(Block block, double siteRateScale);

        void run(Block block, double siteRateScale) {
            changedSiteStamp++;
            ambiguityLogProbability = sampleAmbiguitySites(block, siteRateScale);
            uniformLogProbability = sampleUniformSites(block, siteRateScale);
        }

        void addChangedSite(int p) {
            if (changedSiteStamps[p] != changedSiteStamp) {
                changedSiteStamps[p] = changedSiteStamp;
                changedSites.add(p);
            }
        }

        /**
         * Scores (or draws) the block states of the variable sites with no events around the block: each
         * keeps its state a throughout the block with probability P(unchanged | a). Draws use thinning, so
         * only the sites that change are visited (in expectation, as well as O(events)).
         *
         * @return the log probability of their (drawn or existing) block states under the full conditional
         */
        private double sampleUniformSites(Block block, double siteRateScale) {
            final double w = delegate.getPatternWeights()[variableSites[0]];
            final double nu = delegate.getSiteRate(variableSites[0], siteRateScale);
            final UniformConditional[] uniform = new UniformConditional[stateCount];
            final double[] probabilityChanged = new double[stateCount];
            double maxProbabilityChanged = 0.0;
            for (int a = 0; a < stateCount; a++) {
                if (uniformCounts[a] > 0) {
                    uniform[a] = new UniformConditional(new Conditional(block, a, w, nu, true), a);
                    probabilityChanged[a] = uniform[a].totalChangeMass;
                    maxProbabilityChanged = Math.max(maxProbabilityChanged, probabilityChanged[a]);
                }
            }

            changedUniformCount = 0;
            final int[] changedCounts = new int[stateCount];
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
                    int a = block.store.getState(referenceNode, p);
                    if (a < 0 || MathUtils.nextDouble() * maxProbabilityChanged >= probabilityChanged[a]) {
                        continue;
                    }
                    logProbability += sampleChangedUniformPattern(block, uniform[a], p);
                    changedCounts[a]++;
                    changedUniformCount++;
                    addChangedSite(p);
                }
            }
            for (int a = 0; a < stateCount; a++) {
                if (uniformCounts[a] > 0) {
                    logProbability += (uniformCounts[a] - changedCounts[a]) * uniform[a].logProbabilityUnchanged;
                }
            }
            return logProbability;
        }
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
        double sampleAmbiguitySites(Block block, double siteRateScale) {
            if (ambiguitySites.length == 0) {
                return 0.0;
            }
            final DefaultMutationStore store = block.store;
            final int blockSize = block.nodes.length;
            final byte[][] blockStates = block.states;

            // the block's neighbours, in the order the conditionals use them
            int neighbourCount = block.topParent == null ? 0 : 1;
            for (int i = 0; i < blockSize; i++) {
                for (int c : block.children[i]) {
                    if (!inBlock[c]) {
                        neighbourCount++;
                    }
                }
            }
            final int[] neighbours = new int[neighbourCount];
            int k = 0;
            if (block.topParent != null) {
                neighbours[k++] = block.topParent.getNumber();
            }
            for (int i = 0; i < blockSize; i++) {
                for (int c : block.children[i]) {
                    if (!inBlock[c]) {
                        neighbours[k++] = c;
                    }
                }
            }
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
            // the groups of the current presence, by state
            final AmbiguityGroup[] groupsByState = new AmbiguityGroup[stateCount];
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
                AmbiguityGroup current = groupsByState[a];
                if (current == null) {
                    probe.set(present, a);
                    current = groups.get(probe);
                    if (current == null) {
                        AmbiguityKey key = new AmbiguityKey();
                        key.set(present.clone(), a);
                        current = new AmbiguityGroup(a, p);
                        groups.put(key, current);
                    }
                    groupsByState[a] = current;
                }
                current.sites.add(p);
            }

            final double w = delegate.getPatternWeights()[variableSites[0]];
            final double nu = delegate.getSiteRate(variableSites[0], siteRateScale);
            for (AmbiguityGroup group : groups.values()) {
                final int count = group.sites.count;
                final Conditional conditional = new Conditional(block, group.representative, w, nu);
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
                // the sites that change, by thinning (here every site in the group has the same probability)
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
                        int p = group.sites.values[(int) index];
                        logProbability += sampleChangedUniformPattern(block, uniform, p);
                        addChangedSite(p);
                        changedCount++;
                    }
                }
                logProbability += (count - changedCount) * uniform.logProbabilityUnchanged;
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
                addChangedSite(p);
            }
        }
    }

    /**
     * The ambiguity-only sites of a path move that share a state and a pattern of ambiguous neighbours (and
     * so a conditional).
     */
    private static final class AmbiguityGroup {
        final int state;
        final int representative;
        final IntBuffer sites = new IntBuffer();

        AmbiguityGroup(int state, int representative) {
            this.state = state;
            this.representative = representative;
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
        final int[] substitutionMarks; // the substitution sites (handled individually) are marked here
        final int substitutionStamp;
        final int hotStamp;            // the hot sites (handled as chains) are marked in hotMarks
        byte[][] neighbourStates;      // P's neighbours' states (set when the ambiguity-only sites are swept)

        LocalPass(int referenceNode, int[] uniformCounts, int[] substitutionMarks, int substitutionStamp,
                  int hotStamp) {
            super(referenceNode, uniformCounts);
            this.substitutionMarks = substitutionMarks;
            this.substitutionStamp = substitutionStamp;
            this.hotStamp = hotStamp;
        }

        /**
         * @return true for a substitution site, a hot site, or one at which one of P's neighbours is ambiguous
         */
        @Override
        boolean isCandidate(int p) {
            if (substitutionMarks[p] == substitutionStamp || hotMarks[p] == hotStamp) {
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
            return substitutionMarks[p] == substitutionStamp || hotMarks[p] == hotStamp ||
                    delegate.isConstantPattern(p);
        }

        @Override
        double sampleAmbiguitySites(Block block, double siteRateScale) {
            final DefaultMutationStore store = block.store;
            final int[] children = block.children[0];

            // P's neighbours, in the order of the path pass
            final int neighbourCount = children.length + (block.topParent == null ? 0 : 1);
            final int[] neighbours = new int[neighbourCount];
            int k = 0;
            if (block.topParent != null) {
                neighbours[k++] = block.topParent.getNumber();
            }
            int childBits = 0;
            int xBit = 0;
            for (int c : children) {
                childBits |= 1 << k;
                if (c == referenceNode) {
                    xBit = 1 << k;
                }
                neighbours[k++] = c;
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
            final LocalAmbiguityGroup[] groupsByState = new LocalAmbiguityGroup[stateCount];
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
                        Arrays.fill(groupsByState, null);
                        byte[] states = neighbourStates[first];
                        for (int p = start; p < end; p++) {
                            if (isExcluded(p)) {
                                continue;
                            }
                            int a = states[p];
                            LocalAmbiguityGroup group = groupsByState[a];
                            if (group == null) {
                                probe.set(present, a);
                                group = groups.get(probe);
                                if (group == null) {
                                    AmbiguityKey key = new AmbiguityKey();
                                    key.set(present.clone(), a);
                                    group = new LocalAmbiguityGroup(a, p, first);
                                    groups.put(key, group);
                                }
                                groupsByState[a] = group;
                                group.addSegment(start, end);
                            }
                            group.count++;
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
                                        addChangedSite(p);
                                    }
                                }
                            } else {
                                for (int p = Math.max(runStart, position); p < runEnd; p++) {
                                    if (!isExcluded(p)) {
                                        block.setState(0, p, neighbourStates[first][p]);
                                        addChangedSite(p);
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

            final double w = delegate.getPatternWeights()[variableSites[0]];
            final double nu = delegate.getSiteRate(variableSites[0], siteRateScale);
            double logProbability = 0.0;
            for (Map.Entry<AmbiguityKey, LocalAmbiguityGroup> entry : groups.entrySet()) {
                final LocalAmbiguityGroup group = entry.getValue();
                final int groupMask = (int) entry.getKey().present[0];
                // these sites are not uniform: take them out of the counts of the uniform sites (read at X)
                if ((groupMask & xBit) != 0) {
                    uniformCounts[group.state] -= group.count;
                }
                final Conditional conditional = new Conditional(block, group.representative, w, nu);
                if (!conditional.present[0]) {
                    // no data below P: its states were cleared above
                    continue;
                }
                final UniformConditional uniform = new UniformConditional(conditional, group.state);
                if (!block.sample) {
                    logProbability += group.count * uniform.logProbabilityUnchanged;
                    continue;
                }
                // the sites that change, by thinning (every site in the group has the same probability)
                final double probabilityChanged = uniform.totalChangeMass;
                int changedCount = 0;
                if (probabilityChanged > 0.0) {
                    final double logSkip = Math.log1p(-probabilityChanged);
                    long index = -1;
                    group.resetCursor();
                    while (true) {
                        index = nextThinned(index, logSkip);
                        if (index >= group.count) {
                            break;
                        }
                        int p = locate(group, (int) index);
                        logProbability += sampleChangedUniformPattern(block, uniform, p);
                        addChangedSite(p);
                        changedCount++;
                    }
                }
                logProbability += (group.count - changedCount) * uniform.logProbabilityUnchanged;
            }
            return logProbability;
        }

        /**
         * @return the site of the group with the given index (in pattern order); indices must be asked for in
         * increasing order after resetCursor
         */
        private int locate(LocalAmbiguityGroup group, int index) {
            final byte[] states = neighbourStates[group.firstPresent];
            while (true) {
                int segmentStart = group.segments.values[2 * group.cursorSegment];
                int segmentEnd = group.segments.values[2 * group.cursorSegment + 1];
                if (group.cursorSite < segmentStart) {
                    group.cursorSite = segmentStart;
                }
                for (; group.cursorSite < segmentEnd; group.cursorSite++) {
                    int p = group.cursorSite;
                    if (states[p] == group.state && !isExcluded(p)) {
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
        final int representative;
        final int firstPresent;   // the neighbour the state is read from
        final IntBuffer segments = new IntBuffer(); // [start, end) of each segment with sites in the group
        int count = 0;
        int cursorSegment;
        int cursorSite;
        int cursorIndex;

        LocalAmbiguityGroup(int state, int representative, int firstPresent) {
            this.state = state;
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

    // ---------------------------------------------------------------------------------------------
    // Conditionals
    // ---------------------------------------------------------------------------------------------

    /**
     * The full conditional distribution of a block's states for one pattern, factorised as the
     * distribution of the top node and of each other node given its parent. The upward (filtering)
     * pass is done on construction; the distribution of a node given its parent's state is computed
     * when first needed.
     * <p>
     * Without per-site rates the filtering is done in linear space and the distributions are
     * normalised in linear space; their logarithms are only taken when asked for.
     */
    private final class Conditional {
        final Block block;
        final double w;
        final double nu;
        // if not negative, the state of every node around the block (for sites with no events there)
        private final int uniformState;
        final boolean[] present;
        private double[][] partials;    // density of everything below block node i given its state, scaled (linear filtering)
        private double[][] logPartials; // log density of everything below block node i given its state (log filtering)
        private double[][] linearMatrices;
        final int[] branchOf;           // the branch above each block node
        private double[] topLogP;
        private double[] topP;
        private final double[][][] logRows;
        private final double[][][] rows;

        Conditional(Block block, int p, double w, double nu) {
            this(block, p, -1, w, nu);
        }

        /**
         * The conditional for a site with state a at every node around the block.
         */
        Conditional(Block block, int a, double w, double nu, boolean uniform) {
            this(block, -1, a, w, nu);
        }

        private Conditional(Block block, int p, int uniformState, double w, double nu) {
            final int[][] children = block.children;
            final int blockSize = children.length;
            this.block = block;
            this.w = w;
            this.nu = nu;
            this.uniformState = uniformState;
            this.present = new boolean[blockSize];
            this.branchOf = new int[blockSize];
            this.logRows = new double[blockSize][][];
            this.rows = new double[blockSize][][];
            for (int i = 0; i < blockSize; i++) {
                for (int j = 0; j < children[i].length; j++) {
                    int c = children[i][j];
                    if (inBlock[c]) {
                        branchOf[block.index[c]] = block.childBranch[i][j];
                    }
                }
            }
            branchOf[blockSize - 1] = block.topBranch;

            block.ensureLogF(nu);
            if (delegate.hasSiteRates() || !filterLinear(p)) {
                partials = null;
                filterLog(p);
            }

            int topIndex = blockSize - 1;
            if (present[topIndex]) {
                int parentState = -1;
                if (block.topParent != null) {
                    parentState = stateOf(block.topParent.getNumber(), p);
                    if (parentState < 0) {
                        throw new IllegalStateException("Parent of the resampling block is ambiguous where the block is not");
                    }
                }
                if (partials != null) {
                    double[] weights = new double[stateCount];
                    double[] partial = partials[topIndex];
                    if (parentState < 0) {
                        for (int s = 0; s < stateCount; s++) {
                            weights[s] = partial[s] * Math.exp(w * block.logFrequencies[s]);
                        }
                    } else {
                        double[] matrix = linearMatrices[block.topBranch];
                        for (int s = 0; s < stateCount; s++) {
                            weights[s] = partial[s] * matrix[parentState * stateCount + s];
                        }
                    }
                    topP = normaliseLinear(weights);
                } else {
                    double[] logWeights = new double[stateCount];
                    if (parentState < 0) {
                        for (int s = 0; s < stateCount; s++) {
                            logWeights[s] = logPartials[topIndex][s] + w * block.logFrequencies[s];
                        }
                    } else {
                        double[] f = block.logF[block.topBranch];
                        for (int s = 0; s < stateCount; s++) {
                            logWeights[s] = logPartials[topIndex][s] + w * f[parentState * stateCount + s];
                        }
                    }
                    topLogP = normaliseLog(logWeights);
                    topP = exp(topLogP);
                }
            }
        }

        private int stateOf(int node, int p) {
            return uniformState >= 0 ? uniformState : block.store.getState(node, p);
        }

        /**
         * Filtering in linear space with exp(w logF) matrices (no exponentials or logarithms per pattern).
         *
         * @return false if the partials underflow, in which case filterLog must be used
         */
        private boolean filterLinear(int p) {
            final double[][] m = block.getLinearMatrices(p >= 0 ? weightClass[p] : weightClass[variableSites[0]]);
            final int[][] children = block.children;
            final int[][] childBranch = block.childBranch;
            final int[] index = block.index;
            final boolean[] inBlock = MutationHistoryResampler.this.inBlock;
            final int blockSize = children.length;
            partials = new double[blockSize][stateCount];
            linearMatrices = m;
            for (int i = 0; i < blockSize; i++) {
                double[] partial = partials[i];
                Arrays.fill(partial, 1.0);
                for (int j = 0; j < children[i].length; j++) {
                    int c = children[i][j];
                    double[] matrix = m[childBranch[i][j]];
                    if (inBlock[c]) {
                        int ci = index[c];
                        if (present[ci]) {
                            present[i] = true;
                            double[] childPartial = partials[ci];
                            for (int s = 0; s < stateCount; s++) {
                                double sum = 0.0;
                                for (int t = 0; t < stateCount; t++) {
                                    sum += matrix[s * stateCount + t] * childPartial[t];
                                }
                                partial[s] *= sum;
                            }
                        }
                    } else {
                        int t = stateOf(c, p);
                        if (t >= 0) {
                            present[i] = true;
                            for (int s = 0; s < stateCount; s++) {
                                partial[s] *= matrix[s * stateCount + t];
                            }
                        }
                    }
                }
                if (present[i]) {
                    double max = 0.0;
                    for (int s = 0; s < stateCount; s++) {
                        max = Math.max(max, partial[s]);
                    }
                    if (!(max > 0.0) || Double.isInfinite(max)) {
                        Arrays.fill(present, false);
                        return false;
                    }
                    double scale = 1.0 / max;
                    for (int s = 0; s < stateCount; s++) {
                        partial[s] *= scale;
                    }
                }
            }
            return true;
        }

        /**
         * Filtering in log space (used with per-site rates or if linear filtering underflows).
         */
        private void filterLog(int p) {
            final int[][] children = block.children;
            final int[][] childBranch = block.childBranch;
            final int[] index = block.index;
            final double[][] logF = block.logF;
            final boolean[] inBlock = MutationHistoryResampler.this.inBlock;
            final int blockSize = children.length;
            final double[] terms = new double[stateCount];
            logPartials = new double[blockSize][stateCount];
            for (int i = 0; i < blockSize; i++) {
                double[] partial = logPartials[i];
                for (int j = 0; j < children[i].length; j++) {
                    int c = children[i][j];
                    double[] f = logF[childBranch[i][j]];
                    if (inBlock[c]) {
                        int ci = index[c];
                        if (present[ci]) {
                            present[i] = true;
                            for (int s = 0; s < stateCount; s++) {
                                for (int t = 0; t < stateCount; t++) {
                                    terms[t] = w * f[s * stateCount + t] + logPartials[ci][t];
                                }
                                partial[s] += logSumExp(terms);
                            }
                        }
                    } else {
                        int t = stateOf(c, p);
                        if (t >= 0) {
                            present[i] = true;
                            for (int s = 0; s < stateCount; s++) {
                                partial[s] += w * f[s * stateCount + t];
                            }
                        }
                    }
                }
            }
        }

        /**
         * @return the distribution of (non-top) block node i given its parent has state s
         */
        double[] row(int i, int s) {
            if (rows[i] == null) {
                rows[i] = new double[stateCount][];
                logRows[i] = new double[stateCount][];
            }
            if (rows[i][s] == null) {
                if (partials != null) {
                    double[] matrix = linearMatrices[branchOf[i]];
                    double[] partial = partials[i];
                    double[] weights = new double[stateCount];
                    for (int t = 0; t < stateCount; t++) {
                        weights[t] = matrix[s * stateCount + t] * partial[t];
                    }
                    rows[i][s] = normaliseLinear(weights);
                } else {
                    block.ensureLogF(nu);
                    double[] f = block.logF[branchOf[i]];
                    double[] terms = new double[stateCount];
                    for (int t = 0; t < stateCount; t++) {
                        terms[t] = w * f[s * stateCount + t] + logPartials[i][t];
                    }
                    logRows[i][s] = normaliseLog(terms);
                    rows[i][s] = exp(logRows[i][s]);
                }
            }
            return rows[i][s];
        }

        /**
         * @return the log distribution of (non-top) block node i given its parent has state s
         */
        double[] logRow(int i, int s) {
            double[] row = row(i, s);
            if (logRows[i][s] == null) {
                logRows[i][s] = log(row);
            }
            return logRows[i][s];
        }

        double[] topP() {
            return topP;
        }

        double[] topLogP() {
            if (topLogP == null) {
                topLogP = log(topP);
            }
            return topLogP;
        }

        double[] logRowOrTop(int i, int s) {
            return i == present.length - 1 ? topLogP() : logRow(i, s);
        }

        double[] rowOrTop(int i, int s) {
            return i == present.length - 1 ? topP : row(i, s);
        }
    }

    /**
     * For a pattern whose unambiguous nodes all have state a throughout the block's neighbourhood: the
     * probability that the block's unambiguous nodes all keep a, and for each node (in pre-order) the
     * probability that it is the first to differ from a. Computed in linear space with one logarithm.
     */
    private static final class UniformConditional {
        final Conditional conditional;
        final int state;
        final double logProbabilityUnchanged;
        final double probabilityUnchanged;
        final double[] firstChangeMass; // probability that this pre-order position is the first to differ from a
        final double totalChangeMass;   // the probability that any node differs from a

        UniformConditional(Conditional conditional, int a) {
            this.conditional = conditional;
            this.state = a;
            int blockSize = conditional.present.length;
            firstChangeMass = new double[blockSize];
            double product = 1.0;
            double total = 0.0;
            for (int position = 0; position < blockSize; position++) {
                int i = blockSize - 1 - position;
                if (!conditional.present[i]) {
                    // ambiguous nodes have no state to change
                    continue;
                }
                double[] row = conditional.rowOrTop(i, a);
                double changed = 0.0;
                for (int t = 0; t < row.length; t++) {
                    if (t != a) {
                        changed += row[t];
                    }
                }
                firstChangeMass[position] = product * changed;
                total += firstChangeMass[position];
                product *= row[a];
            }
            probabilityUnchanged = product;
            totalChangeMass = total;
            if (total < 0.5) {
                // more accurate than log(product) when changes are rare (the usual case)
                logProbabilityUnchanged = Math.log1p(-total);
            } else if (product > 1E-250 || Double.isNaN(product)) {
                logProbabilityUnchanged = Math.log(product);
            } else {
                // the product may have underflowed: sum the logarithms
                double logProduct = 0.0;
                for (int position = 0; position < blockSize; position++) {
                    int i = blockSize - 1 - position;
                    if (conditional.present[i]) {
                        logProduct += conditional.logRowOrTop(i, a)[a];
                    }
                }
                logProbabilityUnchanged = logProduct;
            }
        }
    }

    /**
     * Draws the block's states for a uniform pattern conditional on at least one of them differing
     * from the surrounding state: picks the first node (in pre-order) to differ, draws its new state,
     * then draws the nodes after it from their usual conditionals. The block's states must already be
     * in the unchanged configuration.
     *
     * @return the log probability of the drawn states under the full conditional
     */
    private double sampleChangedUniformPattern(Block block, UniformConditional uniform, int p) {
        final Conditional conditional = uniform.conditional;
        final int blockSize = conditional.present.length;
        final int a = uniform.state;

        double target = MathUtils.nextDouble() * uniform.totalChangeMass;
        int first = -1;
        for (int position = 0; position < blockSize; position++) {
            if (uniform.firstChangeMass[position] > 0.0) {
                first = position;
                if (target < uniform.firstChangeMass[position]) {
                    break;
                }
                target -= uniform.firstChangeMass[position];
            }
        }
        if (first < 0) {
            return Double.NEGATIVE_INFINITY;
        }

        // the nodes before the first to change keep a
        double logProbability = 0.0;
        for (int position = 0; position < first; position++) {
            int j = blockSize - 1 - position;
            if (conditional.present[j]) {
                logProbability += conditional.logRowOrTop(j, a)[a];
            }
        }

        // the first node to change: a state other than a, in proportion to its conditional
        int i = blockSize - 1 - first;
        double[] logRow = conditional.logRowOrTop(i, a);
        double[] row = conditional.rowOrTop(i, a);
        double changed = 0.0;
        for (int t = 0; t < stateCount; t++) {
            if (t != a) {
                changed += row[t];
            }
        }
        double u = MathUtils.nextDouble() * changed;
        int newState = -1;
        for (int t = 0; t < stateCount; t++) {
            if (t != a && row[t] > 0.0) {
                newState = t;
                if (u < row[t]) {
                    break;
                }
                u -= row[t];
            }
        }
        if (newState < 0) {
            return Double.NEGATIVE_INFINITY;
        }
        block.setState(i, p, (byte) newState);
        logProbability += logRow[newState];

        // later nodes in pre-order (whose parents are all already drawn) from their usual conditionals
        for (int position = first + 1; position < blockSize; position++) {
            int j = blockSize - 1 - position;
            if (!conditional.present[j]) {
                continue;
            }
            int parentState = block.states[block.parentIndex[j]][p];
            logProbability += choose(conditional.logRow(j, parentState), conditional.row(j, parentState), true, -1);
            block.setState(j, p, (byte) chosenState);
        }
        return logProbability;
    }

    private static final class ConditionalKey {
        final byte[] boundary;
        final double weight;
        final double rate;
        final int hash;

        ConditionalKey(byte[] boundary, double weight, double rate) {
            this.boundary = boundary;
            this.weight = weight;
            this.rate = rate;
            this.hash = 31 * (31 * Arrays.hashCode(boundary) + Double.hashCode(weight)) + Double.hashCode(rate);
        }

        @Override
        public int hashCode() {
            return hash;
        }

        @Override
        public boolean equals(Object o) {
            if (!(o instanceof ConditionalKey)) {
                return false;
            }
            ConditionalKey other = (ConditionalKey) o;
            return weight == other.weight && rate == other.rate && Arrays.equals(boundary, other.boundary);
        }
    }

    /**
     * Draws a state from a normalised distribution or, if not sampling, scores the existing state.
     *
     * @return the log probability of the chosen state, which is left in chosenState
     */
    private double choose(double[] logP, double[] p, boolean sample, int state) {
        if (!sample) {
            if (state < 0) {
                throw new IllegalStateException("MAT has an ambiguous state at a node that has data below it");
            }
            chosenState = state;
            return Double.isNaN(logP[state]) ? Double.NEGATIVE_INFINITY : logP[state];
        }
        if (Double.isNaN(logP[0])) {
            // the conditional has no mass: the proposal is impossible
            chosenState = 0;
            return Double.NEGATIVE_INFINITY;
        }
        double u = MathUtils.nextDouble();
        double cumulative = 0.0;
        for (int s = 0; s < p.length; s++) {
            cumulative += p[s];
            if (u < cumulative) {
                chosenState = s;
                return logP[s];
            }
        }
        // rounding: take the last state with non-zero probability
        for (int s = p.length - 1; s >= 0; s--) {
            if (p[s] > 0.0) {
                chosenState = s;
                return logP[s];
            }
        }
        chosenState = 0;
        return Double.NEGATIVE_INFINITY;
    }

    private static double[] normaliseLog(double[] logWeights) {
        double total = logSumExp(logWeights);
        double[] result = new double[logWeights.length];
        for (int s = 0; s < logWeights.length; s++) {
            result[s] = logWeights[s] - total;
        }
        return result;
    }

    /**
     * Normalises in place; with no mass the result is NaN throughout (as normaliseLog gives).
     */
    private static double[] normaliseLinear(double[] weights) {
        double total = 0.0;
        for (double weight : weights) {
            total += weight;
        }
        if (!(total > 0.0)) {
            Arrays.fill(weights, Double.NaN);
            return weights;
        }
        double scale = 1.0 / total;
        for (int s = 0; s < weights.length; s++) {
            weights[s] *= scale;
        }
        return weights;
    }

    private static double[] exp(double[] logP) {
        double[] p = new double[logP.length];
        for (int s = 0; s < logP.length; s++) {
            p[s] = Math.exp(logP[s]);
        }
        return p;
    }

    private static double[] log(double[] p) {
        double[] logP = new double[p.length];
        for (int s = 0; s < p.length; s++) {
            logP[s] = Math.log(p[s]);
        }
        return logP;
    }

    private static double logSumExp(double[] x) {
        double max = Double.NEGATIVE_INFINITY;
        for (double v : x) {
            if (v > max) {
                max = v;
            }
        }
        if (max == Double.NEGATIVE_INFINITY) {
            return max;
        }
        double sum = 0.0;
        for (double v : x) {
            sum += Math.exp(v - max);
        }
        return max + Math.log(sum);
    }
}
