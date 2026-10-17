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
 * Resamples the mutational history of a MAT jointly with a subtree prune-and-regraft move.
 *
 * With at most one mutation per site per branch and mutation times integrated out, a mutation
 * history is just the state of every pattern at every internal node. When the subtree below
 * {@code node} moves from beside {@code sibling} to the branch above {@code destination}, the only
 * states that need to change are those of the regrafted parent node and of the internal nodes on the
 * path between {@code sibling} and {@code destination} in the pruned tree. This block is the same for
 * the move and its reverse, and contains every node whose ambiguity status can change.
 *
 * The block's states are drawn from their exact full conditional given all other states (a Gibbs
 * block update) by forward filtering and backward sampling on the block, which is a path. The
 * proposal density of the old states (under the old topology) and of the new states (under the new
 * one) gives the Hastings ratio, so the combined move is a valid Metropolis-Hastings move whatever
 * topology move it accompanies.
 *
 * Only sites that vary around the block are resampled (a hybrid of Gibbs sampling and keeping the
 * states). A site is resampled if the states of the block's neighbours (the nodes adjacent to the block
 * in either the old or the new topology) are not all the same; every other site keeps its states, a
 * proposal that is its own reverse. The move never changes the states outside the block, so the move
 * and its reverse resample the same sites. Candidate sites are found from the events on the branches
 * around the block, so a move costs O(events around the block), not O(sites).
 *
 * The hybrid cannot create back-mutations at sites that look uniform around every block it sees, so
 * the alternative {@link SiteSelection#ALL} resamples every (non-constant) site, a plain Gibbs block
 * update. Without per-site rates or explicit hits, it is done sparsely (as in Delphy): a site with no
 * events around the block has the same state a throughout, and its block states stay unchanged with a
 * probability that depends only on a. So the sites that change are drawn by thinning (stepping
 * through the variable sites with geometric jumps at the highest probability of a change, and keeping
 * each site visited with its own probability over that), the unchanged ones contribute
 * n_a log P(unchanged | a) for the number n_a of them in each state, and a move costs
 * O(events around the block + changes), not O(sites). Sites whose only events around the block are
 * ambiguities are treated in the same way, grouped by their state and by which of the block's
 * neighbours are ambiguous there (which together determine their conditional): there are typically
 * many such sites (runs of N) but few groups. Branch history models with explicit hits ({@link BranchHistoryModel#hasExplicitHistories()})
 * require it: the ancestral states are drawn with the hits summed out, then each site's hits on every
 * branch touching the block are drawn given its end states, and the probabilities of the old and new
 * hits enter the Hastings ratio.
 */
public class MutationHistoryResampler implements SubtreeMoveHook {

    /**
     * Which sites a move resamples: those that vary around the block (HYBRID) or all of them (ALL).
     */
    public enum SiteSelection {
        HYBRID, ALL;

        public static SiteSelection parse(String name) {
            for (SiteSelection selection : values()) {
                if (selection.name().equalsIgnoreCase(name)) {
                    return selection;
                }
            }
            throw new IllegalArgumentException("Unknown site selection: " + name);
        }
    }

    private final MATDataLikelihoodDelegate delegate;
    private final BranchHistoryModel historyModel;
    private final SiteSelection siteSelection;
    private double logProposalOfOldHistories;
    private final Tree tree;
    private final int stateCount;
    private final int patternCount;

    /** Pattern classes relative to a block's neighbourhood (see {@link #getPatternClasses}) */
    public static final byte HAS_SUBSTITUTIONS = 0;
    public static final byte NO_EVENTS = 1;
    public static final byte AMBIGUITIES_ONLY = 2;
    public static final byte CONSTANT = 3;

    private int[] block = null;
    private int[] movedNodes = null; // the nodes whose parents a subtree move changes
    private int[] candidateSites = null;
    private int[] oldBoundary = null;

    // per-site work arrays: classes and old/new proposal log probabilities of candidate sites
    private final byte[] siteClasses;
    private final int[] siteStamps;
    private int stamp = 0;
    private final double[] oldSiteLogProbabilities;
    private final double[] newSiteLogProbabilities;
    private final boolean[] resampledSites;
    private int[] lastResampledSites = new int[0];

    private final double[] infinitesimalMatrix;
    private final double[] exitRates;
    private final double[] logFrequencies;
    private final boolean[] inBlock;
    private int chosenState;

    // patterns grouped by weight, so conditionals can be looked up by [state][weight class]
    private final int[] weightClass;
    private final double[] classWeights;

    // sparse handling of the sites with no events around the block (site selection ALL; see the class comment)
    private final boolean sparseUniform;
    private final int[] variableSites;
    private int[][] stateCountsCache;
    private long[] stateCountsVersion;
    private int candidateStamp;
    private int referenceNode;
    private final int[] uniformCounts;
    private double uniformOldLogProbability;
    private int[] changedUniformSites = new int[16];
    private int changedUniformCount;
    // the sites with only ambiguities around the block, handled in groups (see the class comment)
    private int[] substitutionSites;
    private int[] ambiguitySites;
    private double ambiguityOldLogProbability;
    private double lastAmbiguityLogProbability;
    private int[] changedAmbiguitySites = new int[16];
    private int changedAmbiguityCount;
    private int[] changedAmbiguityStamps;
    private int changedAmbiguityStamp = 0;
    // the sites that sampleUniformSites skips (handled individually or in groups)
    private int[] uniformSkipStamps;
    private int uniformSkipStamp;

    // the local block (see LocalMove)
    private final boolean localBlock;
    private LocalMove localMove = null;
    private int[] oldMarks;
    private int[] newMarks;
    private int[] hotMarks;
    private int localStamp = 0;

    public MutationHistoryResampler(MATDataLikelihoodDelegate delegate) {
        this(delegate, null);
    }

    public MutationHistoryResampler(MATDataLikelihoodDelegate delegate, SiteSelection siteSelection) {
        this(delegate, siteSelection, false);
    }

    /**
     * @param siteSelection which sites to resample (null for the default: HYBRID unless the branch
     *                      history model has explicit histories)
     * @param localBlock    resample the local block (the regrafted parent, and the nodes whose data at a site
     *                      come only from the moved subtree) instead of the path block
     */
    public MutationHistoryResampler(MATDataLikelihoodDelegate delegate, SiteSelection siteSelection,
                                    boolean localBlock) {
        this.delegate = delegate;
        this.historyModel = delegate.getHistoryModel();
        if (siteSelection == null) {
            siteSelection = historyModel.hasExplicitHistories() ? SiteSelection.ALL : SiteSelection.HYBRID;
        }
        if (siteSelection == SiteSelection.HYBRID && historyModel.hasExplicitHistories()) {
            throw new IllegalArgumentException("The " + historyModel.getName() + " branch history model needs every site " +
                    "to be resampled (site selection " + SiteSelection.ALL + ")");
        }
        this.siteSelection = siteSelection;
        this.tree = delegate.getTree();
        this.stateCount = delegate.getPatternList().getDataType().getStateCount();
        this.patternCount = delegate.getPatternList().getPatternCount();
        this.infinitesimalMatrix = new double[stateCount * stateCount];
        this.exitRates = new double[stateCount];
        this.logFrequencies = new double[stateCount];
        this.inBlock = new boolean[tree.getNodeCount()];
        this.siteClasses = new byte[patternCount];
        this.siteStamps = new int[patternCount];
        this.oldSiteLogProbabilities = new double[patternCount];
        this.newSiteLogProbabilities = new double[patternCount];
        this.resampledSites = new boolean[patternCount];

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
        boolean equalWeights = true;
        for (int p = 0; p < patternCount; p++) {
            if (!delegate.isConstantPattern(p)) {
                if (variableCount > 0 && weights[p] != weights[variable[0]]) {
                    equalWeights = false;
                }
                variable[variableCount++] = p;
            }
        }
        this.variableSites = Arrays.copyOf(variable, variableCount);
        this.uniformCounts = new int[stateCount];
        this.sparseUniform = siteSelection == SiteSelection.ALL && !historyModel.hasExplicitHistories() &&
                !delegate.hasSiteRates() && equalWeights && variableCount > 0;
        this.localBlock = localBlock;
        if (localBlock && !sparseUniform) {
            throw new IllegalArgumentException("The local block needs site selection " + SiteSelection.ALL +
                    ", a branch history model without explicit hits, no per-site rates and equal pattern weights");
        }
    }

    public MATDataLikelihoodDelegate getDelegate() {
        return delegate;
    }

    @Override
    public void beforeSubtreeMove(NodeRef node, NodeRef sibling, NodeRef destination) {
        if (localBlock) {
            beforeLocalMove(node, sibling, destination);
            return;
        }
        delegate.updateMAT();
        delegate.ensureRates();
        movedNodes = new int[]{sibling.getNumber(), tree.getParent(node).getNumber(), destination.getNumber()};
        block = getBlock(tree.getParent(node), sibling, destination);
        if (sparseUniform) {
            // the sites with events around the block; the others are handled together (see the class comment)
            candidateSites = getCandidateSites(block);
            candidateStamp = stamp;
            uniformSkipStamps = siteStamps;
            uniformSkipStamp = stamp;
            referenceNode = node.getNumber();
            countUniformSites();
            int substitutionCount = 0;
            int ambiguityCount = 0;
            int[] substitution = new int[candidateSites.length];
            int[] ambiguity = new int[candidateSites.length];
            for (int p : candidateSites) {
                if (delegate.isConstantPattern(p)) {
                    continue;
                }
                if (siteClasses[p] == AMBIGUITIES_ONLY) {
                    ambiguity[ambiguityCount++] = p;
                } else {
                    substitution[substitutionCount++] = p;
                }
            }
            substitutionSites = Arrays.copyOf(substitution, substitutionCount);
            // in pattern order, so that they can be grouped by sweeping the neighbours' ambiguous runs
            ambiguitySites = sortSites(ambiguity, ambiguityCount);
        } else {
            candidateSites = siteSelection == SiteSelection.ALL ? getVariableSites(block) : getCandidateSites(block);
        }
        oldBoundary = getBoundary(block);
        // score the candidates' current states under the old topology (only the resampled ones will count)
        uniformPass = sparseUniform;
        sampleBlock(block, false, siteClasses, sparseUniform ? substitutionSites : candidateSites, oldSiteLogProbabilities);
        uniformPass = false;
        uniformOldLogProbability = sparseUniform ? lastUniformLogProbability : 0.0;
        ambiguityOldLogProbability = sparseUniform ? lastAmbiguityLogProbability : 0.0;
        if (historyModel.hasExplicitHistories()) {
            // every candidate is resampled (site selection ALL), so score all their hits around the block
            logProposalOfOldHistories = scoreHistories(getIncidentBranches(block), candidateSites);
        }
    }

    @Override
    public double afterSubtreeMove() {
        if (localBlock) {
            return afterLocalMove();
        }
        if (block == null) {
            throw new IllegalStateException("afterSubtreeMove called without beforeSubtreeMove");
        }
        final int[] currentBlock = block;
        block = null;

        final DefaultMutationStore store = delegate.getMutableMutationStore();
        final List<Integer> reparented = new ArrayList<Integer>();
        IntBuffer synced = syncParents(movedNodes);
        for (int i = 0; i < synced.count; i++) {
            reparented.add(synced.values[i]);
        }

        // the sites to resample: candidates whose states on the block's neighbours (old and new) are not all equal
        final int[] newBoundary = getBoundary(currentBlock);
        final int[] sites = new int[candidateSites.length];
        int siteCount = 0;
        for (int p : sparseUniform ? substitutionSites : candidateSites) {
            int state = store.getState(oldBoundary[0], p);
            boolean uniform = true;
            for (int n : oldBoundary) {
                uniform &= store.getState(n, p) == state;
            }
            for (int n : newBoundary) {
                uniform &= store.getState(n, p) == state;
            }
            if (!uniform || siteSelection == SiteSelection.ALL) {
                sites[siteCount++] = p;
                resampledSites[p] = true;
            }
        }
        final int[] resampled = Arrays.copyOf(sites, siteCount);
        lastResampledSites = resampled;

        double logProposalOfOldStates = uniformOldLogProbability + ambiguityOldLogProbability;
        for (int p : resampled) {
            logProposalOfOldStates += oldSiteLogProbabilities[p];
        }
        uniformPass = sparseUniform;
        changedUniformCount = 0;
        changedAmbiguityCount = 0;
        double logProposalOfNewStates = sampleBlock(currentBlock, true, siteClasses, resampled, newSiteLogProbabilities);
        uniformPass = false;
        delegate.setHistoryResampled();

        int[] changed = resampled;
        if (sparseUniform && changedUniformCount + changedAmbiguityCount > 0) {
            // the other sites whose block states changed also need their branches updated
            changed = Arrays.copyOf(resampled, resampled.length + changedAmbiguityCount + changedUniformCount);
            int k = resampled.length;
            for (int i = 0; i < changedAmbiguityCount; i++) {
                changed[k++] = changedAmbiguitySites[i];
                resampledSites[changedAmbiguitySites[i]] = true;
            }
            for (int i = 0; i < changedUniformCount; i++) {
                changed[k++] = changedUniformSites[i];
                resampledSites[changedUniformSites[i]] = true;
            }
            lastResampledSites = changed;
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
            oldBoundary = null;
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
                if (!fullUpdates.contains(c) && !inBlock(currentBlock, c)) {
                    store.updateBranchSites(c, resampledSites, changed);
                }
            }
        }
        for (int p : changed) {
            resampledSites[p] = false;
        }
        candidateSites = null;
        oldBoundary = null;
        delegate.markTopologySynchronised();

        return logProposalOfOldStates - logProposalOfNewStates;
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
     * runs of ambiguous patterns of the neighbours, so a move costs O(events near P + runs + changes),
     * whatever the distance moved (see optimization_notes/mat_local_spr/DESIGN.md).
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
        int[] newAmbiguitySites;       // in pattern order, including any new hot sites
        int[] newUniformCounts;
        int newStamp;
        double logProposalOfOldStates;
    }

    private static final class IntBuffer {
        int[] values = new int[16];
        int count = 0;

        void add(int value) {
            if (count == values.length) {
                values = Arrays.copyOf(values, 2 * count);
            }
            values[count++] = value;
        }

        int[] toArray() {
            return Arrays.copyOf(values, count);
        }
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
        move.newParent = newParentNode == null ? -1 : newParentNode.getNumber();
        final int s = sibling.getNumber();
        final int d = destination.getNumber();
        move.s = s;
        move.d = d;

        move.presentX = complementRuns(getAmbiguousRuns(move.x));

        // The old chain: G and each ancestor above it lose their states where X was their only data, up to
        // the first node that is above X's new position too.
        final int newAncestorStamp = nextNodeStamp();
        for (NodeRef n = newParentNode; n != null; n = tree.getParent(n)) {
            nodeMarks[n.getNumber()] = newAncestorStamp;
        }
        final IntBuffer chain = new IntBuffer();
        final List<int[]> levels = new ArrayList<int[]>();
        int[] hot = intersectRuns(move.presentX, getAmbiguousRuns(s));
        int previous = move.p;
        for (NodeRef n = gNode; n != null && hot.length > 0 && nodeMarks[n.getNumber()] != newAncestorStamp; n = tree.getParent(n)) {
            for (int i = 0; i < tree.getChildCount(n); i++) {
                int c = tree.getChild(n, i).getNumber();
                if (c != previous) {
                    hot = intersectRuns(hot, getAmbiguousRuns(c));
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
        final IntBuffer oldSubstitutions = new IntBuffer();
        final IntBuffer oldAmbiguities = new IntBuffer();
        collectSubstitutionSites(store, gNode == null ? new int[]{move.x, s} : new int[]{move.x, s, move.p},
                oldMarks, markStamp, oldSubstitutions);
        collectAmbiguousSites(gNode == null ? new int[]{move.x, s} : new int[]{gNode.getNumber(), s, move.x},
                oldMarks, markStamp, oldAmbiguities);
        final int[] oldUniformCounts = getUniformCounts(move.x, oldSubstitutions, oldAmbiguities);

        // The new pass: P with neighbours N, D and X. P's state there is affected by the substitutions
        // between X and N (in the old topology, where they are on the path between them) and between D and
        // N, by the runs of ambiguous patterns of X, D and N, and by the old chain (D may be in it).
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
            pathBranches.add(d);
        }
        final IntBuffer newSubstitutions = new IntBuffer();
        final IntBuffer newAmbiguities = new IntBuffer();
        collectSubstitutionSites(store, pathBranches.toArray(), newMarks, markStamp, newSubstitutions);
        collectAmbiguousSites(newParentNode == null ? new int[]{d, move.x} : new int[]{move.newParent, d, move.x},
                newMarks, markStamp, newAmbiguities);
        for (int site : move.oldHotSites) {
            if (newMarks[site] != markStamp) {
                newMarks[site] = markStamp;
                newAmbiguities.add(site);
            }
        }
        move.newSubstitutionSites = newSubstitutions.toArray();
        move.newAmbiguitySites = sortSites(newAmbiguities.values, newAmbiguities.count);
        move.newUniformCounts = getUniformCounts(move.x, newSubstitutions, newAmbiguities);
        move.newStamp = markStamp;

        // score P's current states (and the old chain's) under the old topology
        referenceNode = move.x;
        uniformSkipStamps = oldMarks;
        uniformSkipStamp = markStamp;
        System.arraycopy(oldUniformCounts, 0, uniformCounts, 0, stateCount);
        substitutionSites = withoutMarked(oldSubstitutions.toArray(), hotMarks, markStamp);
        int[] oldAmbiguitySites = withoutMarked(oldAmbiguities.toArray(), hotMarks, markStamp);
        ambiguitySites = sortSites(oldAmbiguitySites, oldAmbiguitySites.length);
        uniformPass = true;
        double logProbability = sampleBlock(new int[]{move.p}, false, null, substitutionSites, null);
        uniformPass = false;
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
            int[] hot = intersectRuns(move.presentX, getAmbiguousRuns(move.newParent));
            NodeRef n = tree.getNode(move.newParent);
            while (hot.length > 0) {
                chain.add(n.getNumber());
                levels.add(hot);
                NodeRef parent = tree.getParent(n);
                if (parent == null) {
                    break;
                }
                hot = intersectRuns(hot, getAmbiguousRuns(parent.getNumber()));
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
        referenceNode = move.x;
        uniformSkipStamps = newMarks;
        uniformSkipStamp = move.newStamp;
        System.arraycopy(move.newUniformCounts, 0, uniformCounts, 0, stateCount);
        substitutionSites = withoutMarked(move.newSubstitutionSites, hotMarks, newHotStamp);
        ambiguitySites = withoutMarked(move.newAmbiguitySites, hotMarks, newHotStamp);
        changedUniformCount = 0;
        changedAmbiguityCount = 0;
        uniformPass = true;
        double logProposalOfNewStates = sampleBlock(new int[]{move.p}, true, null, substitutionSites, null);
        uniformPass = false;
        logProposalOfNewStates += sampleChains(move.p, newChain, newChainLevels, true);
        delegate.setHistoryResampled();
        // the Fitch sets no longer match the topology; parsimony will recompute them if it is ever needed
        store.setFitchSetsStale(true);

        // the sites at which any state changed
        final IntBuffer changed = new IntBuffer();
        addFlagged(changed, substitutionSites, substitutionSites.length);
        addFlagged(changed, changedAmbiguitySites, changedAmbiguityCount);
        addFlagged(changed, changedUniformSites, changedUniformCount);
        addFlagged(changed, move.oldHotSites, move.oldHotSites.length);
        addFlagged(changed, newHotSites, newHotSites.length);
        final int[] changedSites = sortSites(changed.values, changed.count);
        lastResampledSites = changedSites;

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
        substitutionSites = null;
        ambiguitySites = null;
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
            int[] exact = k + 1 < chain.length ? subtractRuns(levels[k], levels[k + 1]) : levels[k];
            int[] sites = variableSitesIn(exact);
            if (sites.length == 0) {
                continue;
            }
            int[] chainBlock = new int[k + 2];
            chainBlock[0] = p;
            System.arraycopy(chain, 0, chainBlock, 1, k + 1);
            logProbability += sampleBlock(chainBlock, sample, null, sites, null);
        }
        return logProbability;
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

    private void addFlagged(IntBuffer buffer, int[] sites, int count) {
        for (int i = 0; i < count; i++) {
            if (!resampledSites[sites[i]]) {
                resampledSites[sites[i]] = true;
                buffer.add(sites[i]);
            }
        }
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
     * Adds the variable sites at which any of the given nodes is ambiguous (and not yet marked).
     */
    private void collectAmbiguousSites(int[] nodes, int[] marks, int markStamp, IntBuffer out) {
        for (int n : nodes) {
            int[] runs = getAmbiguousRuns(n);
            for (int i = 0; i < runs.length; i += 2) {
                for (int site = runs[i]; site < runs[i + 1]; site++) {
                    if (marks[site] != markStamp && !delegate.isConstantPattern(site)) {
                        marks[site] = markStamp;
                        out.add(site);
                    }
                }
            }
        }
    }

    /**
     * @return the number of variable sites in each state at node x, excluding the given sites
     */
    private int[] getUniformCounts(int x, IntBuffer... excluded) {
        final DefaultMutationStore store = delegate.getMutableMutationStore();
        final int[] counts = getStateCounts(x).clone();
        for (IntBuffer sites : excluded) {
            for (int i = 0; i < sites.count; i++) {
                int a = store.getState(x, sites.values[i]);
                if (a >= 0) {
                    counts[a]--;
                }
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
     * @return the intersection of two sets of runs ([start, end) pairs in pattern order)
     */
    static int[] intersectRuns(int[] a, int[] b) {
        int[] result = new int[a.length + b.length];
        int count = 0;
        int i = 0;
        int j = 0;
        while (i < a.length && j < b.length) {
            int start = Math.max(a[i], b[j]);
            int end = Math.min(a[i + 1], b[j + 1]);
            if (start < end) {
                result[count++] = start;
                result[count++] = end;
            }
            if (a[i + 1] < b[j + 1]) {
                i += 2;
            } else {
                j += 2;
            }
        }
        return Arrays.copyOf(result, count);
    }

    /**
     * @return the runs of patterns not in the given runs
     */
    private int[] complementRuns(int[] runs) {
        int[] result = new int[runs.length + 2];
        int count = 0;
        int position = 0;
        for (int i = 0; i < runs.length; i += 2) {
            if (runs[i] > position) {
                result[count++] = position;
                result[count++] = runs[i];
            }
            position = runs[i + 1];
        }
        if (position < patternCount) {
            result[count++] = position;
            result[count++] = patternCount;
        }
        return Arrays.copyOf(result, count);
    }

    private int[] subtractRuns(int[] a, int[] b) {
        return intersectRuns(a, complementRuns(b));
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

    /**
     * @return the sites resampled by the last move (for testing)
     */
    public int[] getLastResampledSites() {
        return lastResampledSites.clone();
    }

    /**
     * @return the number of sites with no events around the block that the last move changed (for testing)
     */
    public int getLastChangedUniformSiteCount() {
        return changedUniformCount;
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

    private static boolean inBlock(int[] blockNodes, int n) {
        for (int b : blockNodes) {
            if (b == n) {
                return true;
            }
        }
        return false;
    }

    /**
     * @return the sites with an event (substitution or ambiguity) on a branch within the block's
     * neighbourhood (the branches above the block nodes and their children), with their classes set in
     * siteClasses. Every other site has the same state throughout the neighbourhood.
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
     * @return the nodes adjacent to the block in the current topology: children of block nodes that
     * are not in the block, and the parent of the top node
     */
    private int[] getBoundary(int[] blockNodes) {
        for (int n : blockNodes) {
            inBlock[n] = true;
        }
        List<Integer> boundary = new ArrayList<Integer>();
        for (int n : blockNodes) {
            NodeRef node = tree.getNode(n);
            for (int i = 0; i < tree.getChildCount(node); i++) {
                int c = tree.getChild(node, i).getNumber();
                if (!inBlock[c]) {
                    boundary.add(c);
                }
            }
            NodeRef parent = tree.getParent(node);
            if (parent != null && !inBlock[parent.getNumber()]) {
                boundary.add(parent.getNumber());
            }
        }
        clearBlock(blockNodes);
        int[] result = new int[boundary.size()];
        for (int i = 0; i < result.length; i++) {
            result[i] = boundary.get(i);
        }
        return result;
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

    /**
     * Draws (or, if sample is false, scores) the states of the block's nodes from their full
     * conditional distribution given all other states in the current tree.
     *
     * @return the log probability of the drawn (or existing) states under that distribution
     */
    public double sampleBlock(int[] blockNodes, boolean sample) {
        return sampleBlock(blockNodes, sample, null, null, null);
    }

    /**
     * As {@link #sampleBlock(int[], boolean)} but with patterns that have no substitutions around the
     * block (see {@link #getPatternClasses}) handled by a fast exact path.
     */
    public double sampleBlock(int[] blockNodes, boolean sample, byte[] patternClasses) {
        return sampleBlock(blockNodes, sample, patternClasses, null, null);
    }

    /**
     * As {@link #sampleBlock(int[], boolean, byte[])} but only for the given sites (all sites if null),
     * optionally recording each site's log probability.
     */
    private double sampleBlock(int[] blockNodes, boolean sample, byte[] patternClasses,
                               int[] sites, double[] siteLogProbabilities) {
        delegate.ensureRates();
        final DefaultMutationStore store = delegate.getMutableMutationStore();
        final int[] postOrder = getBlockPostOrder(blockNodes);
        final int blockSize = postOrder.length;
        final int top = postOrder[blockSize - 1];
        final NodeRef topParent = tree.getParent(tree.getNode(top));

        for (int n : postOrder) {
            inBlock[n] = true;
        }

        // local indices: the block nodes, and the branches below each of them plus the one above the top
        final int[][] children = new int[blockSize][];
        if (nodeBlockIndex == null) {
            nodeBlockIndex = new int[tree.getNodeCount()];
        }
        // (only read for block nodes, so it need not be cleared)
        final int[] blockIndex = nodeBlockIndex;
        for (int i = 0; i < blockSize; i++) {
            blockIndex[postOrder[i]] = i;
        }
        int branchTotal = 1;
        for (int i = 0; i < blockSize; i++) {
            branchTotal += tree.getChildCount(tree.getNode(postOrder[i]));
        }
        final double[] lengths = new double[branchTotal];
        int branchCounter = 0;
        final int[][] childBranch = new int[blockSize][];
        for (int i = 0; i < blockSize; i++) {
            NodeRef node = tree.getNode(postOrder[i]);
            children[i] = new int[tree.getChildCount(node)];
            childBranch[i] = new int[children[i].length];
            for (int j = 0; j < children[i].length; j++) {
                NodeRef child = tree.getChild(node, j);
                children[i][j] = child.getNumber();
                childBranch[i][j] = branchCounter;
                lengths[branchCounter++] = delegate.getBranchLength(child);
            }
        }
        final int[] parentIndex = new int[blockSize];
        parentIndex[blockSize - 1] = -1;
        for (int i = 0; i < blockSize; i++) {
            for (int c : children[i]) {
                if (inBlock[c]) {
                    parentIndex[blockIndex[c]] = i;
                }
            }
        }
        final int topBranch = branchCounter;
        lengths[branchCounter++] = delegate.getBranchLength(tree.getNode(top));
        final int branchCount = branchCounter;

        // model quantities
        delegate.getSubstitutionModel().getInfinitesimalMatrix(infinitesimalMatrix);
        for (int a = 0; a < stateCount; a++) {
            double q = 0.0;
            for (int b = 0; b < stateCount; b++) {
                if (b != a) {
                    q += infinitesimalMatrix[a * stateCount + b];
                }
            }
            exitRates[a] = q;
        }
        double[] frequencies = delegate.getRootFrequencyModel().getFrequencies();
        for (int a = 0; a < stateCount; a++) {
            logFrequencies[a] = Math.log(frequencies[a]);
        }
        final double siteRateScale = delegate.getSiteRateScale();
        final double[] patternWeights = delegate.getPatternWeights();

        // per-call context used by the conditionals
        this.blockChildren = children;
        this.blockChildBranch = childBranch;
        this.blockIndex = blockIndex;
        this.blockParentIndex = parentIndex;
        this.blockTopParent = topParent;
        this.blockTopBranch = topBranch;
        this.blockStore = store;
        this.blockNodes = postOrder;
        this.branchLengths = lengths;
        this.logF = new double[branchCount][stateCount * stateCount];
        this.logFRate = Double.NaN;
        this.linearMatrices = new double[classWeights.length][][];
        ensureLogF(siteRateScale);

        final byte[][] blockStates = new byte[blockSize][];
        for (int i = 0; i < blockSize; i++) {
            // changed only through setBlockState, which keeps the store's node summaries up to date
            blockStates[i] = sample ? store.getStatesForUpdate(postOrder[i]) : store.getStates(postOrder[i]);
        }

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
        final int topIndex = blockSize - 1;

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
                        uniform = new UniformConditional(new Conditional(p, w, nu), a);
                        uniformByClass[a][weightClass[p]] = uniform;
                    }
                } else {
                    ConditionalKey key = new ConditionalKey(new byte[]{(byte) a}, w, nu);
                    uniform = uniformCache.get(key);
                    if (uniform == null) {
                        uniform = new UniformConditional(new Conditional(p, w, nu), a);
                        uniformCache.put(key, uniform);
                    }
                }
                if (!sample || MathUtils.nextDouble() < uniform.probabilityUnchanged) {
                    logProbability += uniform.logProbabilityUnchanged;
                } else {
                    logProbability += sampleChangedUniformPattern(uniform, p, blockStates);
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
                conditional = new Conditional(p, w, nu);
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
                            setBlockState(blockStates, i, p, -1);
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
                    setBlockState(blockStates, i, p, (byte) (conditional.present[i] ? a : -1));
                }
                if (MathUtils.nextDouble() < uniform.probabilityUnchanged) {
                    logProbability += uniform.logProbabilityUnchanged;
                } else {
                    logProbability += sampleChangedUniformPattern(uniform, p, blockStates);
                }
                continue;
            }

            if (!conditional.present[topIndex]) {
                if (sample) {
                    for (int i = 0; i < blockSize; i++) {
                        setBlockState(blockStates, i, p, -1);
                    }
                }
                continue;
            }
            logProbability += choose(conditional.topLogP(), conditional.topP(), sample, blockStates[topIndex][p]);
            if (sample) {
                setBlockState(blockStates, topIndex, p, (byte) chosenState);
            }
            // pre-order: each child given its (already chosen) parent
            for (int i = topIndex - 1; i >= 0; i--) {
                if (!conditional.present[i]) {
                    if (sample) {
                        setBlockState(blockStates, i, p, -1);
                    }
                    continue;
                }
                int s = blockStates[parentIndex[i]][p];
                logProbability += choose(conditional.logRow(i, s), conditional.row(i, s), sample, blockStates[i][p]);
                if (sample) {
                    setBlockState(blockStates, i, p, (byte) chosenState);
                }
            }
        }

        if (siteLogProbabilities != null && lastSite >= 0) {
            siteLogProbabilities[lastSite] = logProbability - lastLogProbability;
        }

        if (uniformPass) {
            lastAmbiguityLogProbability = sampleAmbiguitySites(sample, blockStates, postOrder, children, topParent,
                    store, siteRateScale);
            logProbability += lastAmbiguityLogProbability;
            lastUniformLogProbability = sampleUniformSites(sample, blockStates, siteRateScale);
            logProbability += lastUniformLogProbability;
        }

        for (int n : postOrder) {
            inBlock[n] = false;
        }
        this.blockStore = null;
        return logProbability;
    }

    private double lastUniformLogProbability;
    // set while sampleBlock is scoring or drawing a move's block (when it also handles the sites with no events)
    private boolean uniformPass = false;

    /**
     * The sites with only ambiguities around the block that share a state and a pattern of ambiguous
     * neighbours (and so a conditional).
     */
    private static final class AmbiguityGroup {
        final int state;
        final int representative;
        int[] sites = new int[8];
        int count = 0;

        AmbiguityGroup(int state, int representative) {
            this.state = state;
            this.representative = representative;
        }

        void add(int p) {
            if (count == sites.length) {
                sites = Arrays.copyOf(sites, 2 * count);
            }
            sites[count++] = p;
        }
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
     * Scores (or draws) the block states of the sites whose only events around the block are
     * ambiguities. Every unambiguous node around and in the block has the same state a at such a site,
     * and its block states stay at a (where the block has data below it) with a probability that depends
     * only on a and on which of the block's neighbours are ambiguous, so the sites are handled in groups
     * as the sites with no events are (see sampleUniformSites).
     *
     * @return the log probability of their (drawn or existing) block states under the full conditional
     */
    private double sampleAmbiguitySites(boolean sample, byte[][] blockStates, int[] postOrder, int[][] children,
                                        NodeRef topParent, DefaultMutationStore store, double siteRateScale) {
        changedAmbiguityCount = 0;
        if (changedAmbiguityStamps == null) {
            changedAmbiguityStamps = new int[patternCount];
        }
        changedAmbiguityStamp++;
        if (ambiguitySites.length == 0) {
            return 0.0;
        }
        final int blockSize = postOrder.length;

        // the block's neighbours, in the order the conditionals use them
        int neighbourCount = topParent == null ? 0 : 1;
        for (int i = 0; i < blockSize; i++) {
            for (int c : children[i]) {
                if (!inBlock[c]) {
                    neighbourCount++;
                }
            }
        }
        final int[] neighbours = new int[neighbourCount];
        int k = 0;
        if (topParent != null) {
            neighbours[k++] = topParent.getNumber();
        }
        for (int i = 0; i < blockSize; i++) {
            for (int c : children[i]) {
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
        int eventCount = 0;
        final int[][] runs = new int[neighbourCount][];
        for (int j = 0; j < neighbourCount; j++) {
            runs[j] = getAmbiguousRuns(neighbours[j]);
            eventCount += runs[j].length;
        }
        final long[] events = new long[eventCount];
        k = 0;
        for (int j = 0; j < neighbourCount; j++) {
            for (int position : runs[j]) {
                events[k++] = ((long) position << 32) | j;
            }
        }
        Arrays.sort(events);

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
                setBlockStates(blockStates, p, -1, null, sample);
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
            current.add(p);
        }

        final double w = delegate.getPatternWeights()[variableSites[0]];
        final double nu = delegate.getSiteRate(variableSites[0], siteRateScale);
        for (AmbiguityGroup group : groups.values()) {
            final Conditional conditional = new Conditional(group.representative, w, nu);
            final int topIndex = blockSize - 1;
            if (!conditional.present[topIndex]) {
                // no data below the block
                for (int i = 0; i < group.count; i++) {
                    setBlockStates(blockStates, group.sites[i], -1, null, sample);
                }
                continue;
            }
            final UniformConditional uniform = new UniformConditional(conditional, group.state);
            if (!sample) {
                logProbability += group.count * uniform.logProbabilityUnchanged;
                continue;
            }
            // the unchanged configuration: unambiguous nodes take a (which they may not have had before the move)
            for (int i = 0; i < group.count; i++) {
                setBlockStates(blockStates, group.sites[i], group.state, conditional.present, true);
            }
            // the sites that change, by thinning (here every site in the group has the same probability)
            final double probabilityChanged = uniform.totalChangeMass;
            int changedCount = 0;
            if (probabilityChanged > 0.0) {
                final double logSkip = Math.log1p(-probabilityChanged);
                long index = -1;
                while (true) {
                    double u = 1.0 - MathUtils.nextDouble(); // in (0, 1]
                    index += 1 + (logSkip == Double.NEGATIVE_INFINITY ? 0L : (long) Math.floor(Math.log(u) / logSkip));
                    if (index >= group.count) {
                        break;
                    }
                    int p = group.sites[(int) index];
                    logProbability += sampleChangedUniformPattern(uniform, p, blockStates);
                    addChangedAmbiguitySite(p);
                    changedCount++;
                }
            }
            logProbability += (group.count - changedCount) * uniform.logProbabilityUnchanged;
        }
        return logProbability;
    }

    /**
     * Sets a site's block states to the state (where present, or everywhere if present is null), noting the
     * site as changed if any of them change. Does nothing if not sampling.
     */
    private void setBlockStates(byte[][] blockStates, int p, int state, boolean[] present, boolean sample) {
        if (!sample) {
            return;
        }
        boolean changed = false;
        for (int i = 0; i < blockStates.length; i++) {
            byte value = (byte) (present == null || present[i] ? state : -1);
            if (blockStates[i][p] != value) {
                setBlockState(blockStates, i, p, value);
                changed = true;
            }
        }
        if (changed) {
            addChangedAmbiguitySite(p);
        }
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
        return sorted;
    }

    /**
     * @return the runs of patterns at which the node is ambiguous, as [start, end) pairs in pattern order
     * (cached for the node's current states)
     */
    private int[] getAmbiguousRuns(int n) {
        final DefaultMutationStore store = delegate.getMutableMutationStore();
        if (ambiguousRuns == null) {
            ambiguousRuns = new int[tree.getNodeCount()][];
            ambiguousRunsVersion = new long[tree.getNodeCount()];
        }
        long version = store.getStateVersion(n);
        if (ambiguousRuns[n] == null || ambiguousRunsVersion[n] != version) {
            byte[] states = store.getStates(n);
            int[] runs = new int[8];
            int count = 0;
            int p = 0;
            while (p < patternCount) {
                if (states[p] < 0) {
                    int start = p;
                    while (p < patternCount && states[p] < 0) {
                        p++;
                    }
                    if (count + 2 > runs.length) {
                        runs = Arrays.copyOf(runs, 2 * runs.length);
                    }
                    runs[count++] = start;
                    runs[count++] = p;
                } else {
                    p++;
                }
            }
            ambiguousRuns[n] = Arrays.copyOf(runs, count);
            ambiguousRunsVersion[n] = version;
        }
        return ambiguousRuns[n];
    }

    /**
     * Sets the state of block node i at a pattern, through the store.
     */
    private void setBlockState(byte[][] blockStates, int i, int p, int state) {
        if (blockStates[i][p] != state) {
            blockStore.setState(blockNodes[i], p, state);
        }
    }

    private void addChangedAmbiguitySite(int p) {
        if (changedAmbiguityStamps[p] == changedAmbiguityStamp) {
            return;
        }
        changedAmbiguityStamps[p] = changedAmbiguityStamp;
        if (changedAmbiguityCount == changedAmbiguitySites.length) {
            changedAmbiguitySites = Arrays.copyOf(changedAmbiguitySites, 2 * changedAmbiguityCount);
        }
        changedAmbiguitySites[changedAmbiguityCount++] = p;
    }

    /**
     * Counts the variable sites with no events around the block, by their state (which is the same at
     * every node around the block, so it is read at the reference node: the root of the moving subtree,
     * whose states the move does not change).
     */
    private void countUniformSites() {
        final DefaultMutationStore store = delegate.getMutableMutationStore();
        final int n = referenceNode;
        System.arraycopy(getStateCounts(n), 0, uniformCounts, 0, stateCount);
        for (int p : candidateSites) {
            if (!delegate.isConstantPattern(p)) {
                int a = store.getState(n, p);
                if (a >= 0) {
                    uniformCounts[a]--;
                }
            }
        }
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
     * Scores (or draws) the block states of the variable sites with no events around the block: each
     * keeps its state a throughout the block with probability P(unchanged | a). Draws use thinning, so
     * only the sites that change are visited (in expectation, as well as O(events)).
     *
     * @return the log probability of their (drawn or existing) block states under the full conditional
     */
    private double sampleUniformSites(boolean sample, byte[][] blockStates, double siteRateScale) {
        final DefaultMutationStore store = blockStore;
        final double w = delegate.getPatternWeights()[variableSites[0]];
        final double nu = delegate.getSiteRate(variableSites[0], siteRateScale);
        final UniformConditional[] uniform = new UniformConditional[stateCount];
        final double[] probabilityChanged = new double[stateCount];
        double maxProbabilityChanged = 0.0;
        for (int a = 0; a < stateCount; a++) {
            if (uniformCounts[a] > 0) {
                uniform[a] = new UniformConditional(new Conditional(a, w, nu, true), a);
                probabilityChanged[a] = uniform[a].totalChangeMass;
                maxProbabilityChanged = Math.max(maxProbabilityChanged, probabilityChanged[a]);
            }
        }

        changedUniformCount = 0;
        final int[] changedCounts = new int[stateCount];
        double logProbability = 0.0;
        if (sample && maxProbabilityChanged > 0.0) {
            final double logSkip = Math.log1p(-maxProbabilityChanged);
            long i = -1;
            while (true) {
                // the number of sites skipped before the next visited one is geometric
                double u = 1.0 - MathUtils.nextDouble(); // in (0, 1]
                i += 1 + (logSkip == Double.NEGATIVE_INFINITY ? 0L : (long) Math.floor(Math.log(u) / logSkip));
                if (i >= variableSites.length) {
                    break;
                }
                int p = variableSites[(int) i];
                if (uniformSkipStamps[p] == uniformSkipStamp) {
                    // a site with events: handled individually
                    continue;
                }
                int a = store.getState(referenceNode, p);
                if (a < 0 || MathUtils.nextDouble() * maxProbabilityChanged >= probabilityChanged[a]) {
                    continue;
                }
                logProbability += sampleChangedUniformPattern(uniform[a], p, blockStates);
                changedCounts[a]++;
                if (changedUniformCount == changedUniformSites.length) {
                    changedUniformSites = Arrays.copyOf(changedUniformSites, 2 * changedUniformCount);
                }
                changedUniformSites[changedUniformCount++] = p;
            }
        }
        for (int a = 0; a < stateCount; a++) {
            if (uniformCounts[a] > 0) {
                logProbability += (uniformCounts[a] - changedCounts[a]) * uniform[a].logProbabilityUnchanged;
            }
        }
        return logProbability;
    }

    // per-call context of sampleBlock
    private int[][] blockChildren;
    private int[][] blockChildBranch;
    private int[] blockIndex;
    private int[] blockParentIndex;
    private NodeRef blockTopParent;
    private int blockTopBranch;
    private DefaultMutationStore blockStore;
    private int[] blockNodes;
    private int[] nodeBlockIndex;
    // stamped marks of nodes (a node is marked if its entry equals nodeStamp)
    private int[] nodeMarks;
    private int nodeStamp = 0;
    private long[] siteBits;
    private int[][] ambiguousRuns;
    private long[] ambiguousRunsVersion;
    private double[] branchLengths;
    private double[][] logF;
    private double logFRate;
    private double[][][] linearMatrices;

    /**
     * Makes logF hold the log branch densities for the given relative rate.
     */
    private void ensureLogF(double nu) {
        if (nu != logFRate) {
            for (int e = 0; e < logF.length; e++) {
                computeLogBranchDensities(branchLengths[e], nu, logF[e]);
            }
            logFRate = nu;
        }
    }

    /**
     * @return exp(w logF) for every branch, for the weight class's w (at the common rate, i.e. without
     * per-site rates)
     */
    private double[][] getLinearMatrices(int weightClassIndex) {
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
     * The full conditional distribution of the block's states for one pattern, factorised as the
     * distribution of the top node and of each other node given its parent. The upward (filtering)
     * pass is done on construction; the distribution of a node given its parent's state is computed
     * when first needed.
     * <p>
     * Without per-site rates the filtering is done in linear space and the distributions are
     * normalised in linear space; their logarithms are only taken when asked for.
     */
    private final class Conditional {
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

        Conditional(int p, double w, double nu) {
            this(p, -1, w, nu);
        }

        /**
         * The conditional for a site with state a at every node around the block.
         */
        Conditional(int a, double w, double nu, boolean uniform) {
            this(-1, a, w, nu);
        }

        private Conditional(int p, int uniformState, double w, double nu) {
            final int blockSize = blockChildren.length;
            this.w = w;
            this.nu = nu;
            this.uniformState = uniformState;
            this.present = new boolean[blockSize];
            this.branchOf = new int[blockSize];
            this.logRows = new double[blockSize][][];
            this.rows = new double[blockSize][][];
            for (int i = 0; i < blockSize; i++) {
                for (int j = 0; j < blockChildren[i].length; j++) {
                    int c = blockChildren[i][j];
                    if (inBlock[c]) {
                        branchOf[blockIndex[c]] = blockChildBranch[i][j];
                    }
                }
            }
            branchOf[blockSize - 1] = blockTopBranch;

            ensureLogF(nu);
            if (delegate.hasSiteRates() || !filterLinear(p)) {
                partials = null;
                filterLog(p);
            }

            int topIndex = blockSize - 1;
            if (present[topIndex]) {
                int parentState = -1;
                if (blockTopParent != null) {
                    parentState = stateOf(blockTopParent.getNumber(), p);
                    if (parentState < 0) {
                        throw new IllegalStateException("Parent of the resampling block is ambiguous where the block is not");
                    }
                }
                if (partials != null) {
                    double[] weights = new double[stateCount];
                    double[] partial = partials[topIndex];
                    if (parentState < 0) {
                        for (int s = 0; s < stateCount; s++) {
                            weights[s] = partial[s] * Math.exp(w * logFrequencies[s]);
                        }
                    } else {
                        double[] matrix = linearMatrices[blockTopBranch];
                        for (int s = 0; s < stateCount; s++) {
                            weights[s] = partial[s] * matrix[parentState * stateCount + s];
                        }
                    }
                    topP = normaliseLinear(weights);
                } else {
                    double[] logWeights = new double[stateCount];
                    if (parentState < 0) {
                        for (int s = 0; s < stateCount; s++) {
                            logWeights[s] = logPartials[topIndex][s] + w * logFrequencies[s];
                        }
                    } else {
                        double[] f = logF[blockTopBranch];
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
            return uniformState >= 0 ? uniformState : blockStore.getState(node, p);
        }

        /**
         * Filtering in linear space with exp(w logF) matrices (no exponentials or logarithms per pattern).
         *
         * @return false if the partials underflow, in which case filterLog must be used
         */
        private boolean filterLinear(int p) {
            final double[][] m = getLinearMatrices(p >= 0 ? weightClass[p] : weightClass[variableSites[0]]);
            final int blockSize = blockChildren.length;
            partials = new double[blockSize][stateCount];
            linearMatrices = m;
            for (int i = 0; i < blockSize; i++) {
                double[] partial = partials[i];
                Arrays.fill(partial, 1.0);
                for (int j = 0; j < blockChildren[i].length; j++) {
                    int c = blockChildren[i][j];
                    double[] matrix = m[blockChildBranch[i][j]];
                    if (inBlock[c]) {
                        int ci = blockIndex[c];
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
            final int blockSize = blockChildren.length;
            final double[] terms = new double[stateCount];
            logPartials = new double[blockSize][stateCount];
            for (int i = 0; i < blockSize; i++) {
                double[] partial = logPartials[i];
                for (int j = 0; j < blockChildren[i].length; j++) {
                    int c = blockChildren[i][j];
                    double[] f = logF[blockChildBranch[i][j]];
                    if (inBlock[c]) {
                        int ci = blockIndex[c];
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
                    ensureLogF(nu);
                    double[] f = logF[branchOf[i]];
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

    /**
     * For a pattern whose unambiguous nodes all have state a throughout the block's neighbourhood: the
     * probability that the block's unambiguous nodes all keep a, and for each node (in pre-order) the
     * probability that it is the first to differ from a. Computed in linear space with one logarithm.
     */
    private final class UniformConditional {
        final Conditional conditional;
        final int state;
        final double logProbabilityUnchanged;
        final double probabilityUnchanged;
        final double[] prefix;          // probability that all nodes before this pre-order position keep a
        final double[] firstChangeMass; // probability that this position is the first to differ from a
        final double totalChangeMass;   // the probability that any node differs from a

        UniformConditional(Conditional conditional, int a) {
            this.conditional = conditional;
            this.state = a;
            int blockSize = conditional.present.length;
            prefix = new double[blockSize];
            firstChangeMass = new double[blockSize];
            double product = 1.0;
            double total = 0.0;
            for (int position = 0; position < blockSize; position++) {
                int i = blockSize - 1 - position;
                prefix[position] = product;
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
    private double sampleChangedUniformPattern(UniformConditional uniform, int p, byte[][] blockStates) {
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
        setBlockState(blockStates, i, p, (byte) newState);
        logProbability += logRow[newState];

        // later nodes in pre-order (whose parents are all already drawn) from their usual conditionals
        for (int position = first + 1; position < blockSize; position++) {
            int j = blockSize - 1 - position;
            if (!conditional.present[j]) {
                continue;
            }
            int parentState = blockStates[blockParentIndex[j]][p];
            logProbability += choose(conditional.logRow(j, parentState), conditional.row(j, parentState), true, -1);
            setBlockState(blockStates, j, p, (byte) chosenState);
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
     * Log density of going from state a to state b along a branch of length R (expected substitutions),
     * summed over whatever the branch history model does not sample.
     */
    private void computeLogBranchDensities(double length, double nu, double[] out) {
        historyModel.getLogTransitionDensities(length, nu, out);
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
