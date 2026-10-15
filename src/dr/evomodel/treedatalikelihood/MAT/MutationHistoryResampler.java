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
 */
public class MutationHistoryResampler implements SubtreeMoveHook {

    private final MATDataLikelihoodDelegate delegate;
    private final Tree tree;
    private final int stateCount;
    private final int patternCount;

    /** Pattern classes relative to a block's neighbourhood (see {@link #getPatternClasses}) */
    public static final byte HAS_SUBSTITUTIONS = 0;
    public static final byte NO_EVENTS = 1;
    public static final byte AMBIGUITIES_ONLY = 2;
    public static final byte CONSTANT = 3;

    private int[] block = null;
    private byte[] patternClasses = null;
    private double logProposalOfOldStates;

    private final double[] infinitesimalMatrix;
    private final double[] exitRates;
    private final double[] logFrequencies;
    private final boolean[] inBlock;
    private int chosenState;

    // patterns grouped by weight, so conditionals can be looked up by [state][weight class]
    private final int[] weightClass;
    private final double[] classWeights;

    public MutationHistoryResampler(MATDataLikelihoodDelegate delegate) {
        this.delegate = delegate;
        this.tree = delegate.getTree();
        this.stateCount = delegate.getPatternList().getDataType().getStateCount();
        this.patternCount = delegate.getPatternList().getPatternCount();
        this.infinitesimalMatrix = new double[stateCount * stateCount];
        this.exitRates = new double[stateCount];
        this.logFrequencies = new double[stateCount];
        this.inBlock = new boolean[tree.getNodeCount()];

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
    }

    public MATDataLikelihoodDelegate getDelegate() {
        return delegate;
    }

    @Override
    public void beforeSubtreeMove(NodeRef node, NodeRef sibling, NodeRef destination) {
        delegate.updateMAT();
        block = getBlock(tree.getParent(node), sibling, destination);
        patternClasses = getPatternClasses(block);
        logProposalOfOldStates = sampleBlock(block, false, patternClasses);
    }

    @Override
    public double afterSubtreeMove() {
        if (block == null) {
            throw new IllegalStateException("afterSubtreeMove called without beforeSubtreeMove");
        }
        int[] currentBlock = block;
        block = null;

        DefaultMutationStore store = delegate.getMutableMutationStore();
        int[] parents = MATInitializer.getParents(tree);
        List<Integer> reparented = new ArrayList<Integer>();
        for (int n = 0; n < parents.length; n++) {
            if (store.getParent(n) != parents[n]) {
                reparented.add(n);
            }
        }
        store.setParents(parents);

        double logProposalOfNewStates = sampleBlock(currentBlock, true, patternClasses);
        patternClasses = null;
        delegate.setHistoryResampled();

        // the Fitch sets no longer match the topology; parsimony will recompute them if it is ever needed
        store.setFitchSetsStale(true);
        int[][] children = MATInitializer.getChildren(tree);

        Set<Integer> branches = new HashSet<Integer>(reparented);
        for (int n : currentBlock) {
            branches.add(n);
            for (int c : children[n]) {
                branches.add(c);
            }
        }
        for (int n : branches) {
            store.updateBranch(n);
        }

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
        return sampleBlock(blockNodes, sample, null);
    }

    /**
     * As {@link #sampleBlock(int[], boolean)} but with patterns that have no substitutions around the
     * block (see {@link #getPatternClasses}) handled by a fast exact path.
     */
    public double sampleBlock(int[] blockNodes, boolean sample, byte[] patternClasses) {
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
        final int[] blockIndex = new int[tree.getNodeCount()];
        for (int i = 0; i < blockSize; i++) {
            blockIndex[postOrder[i]] = i;
        }
        final List<Double> lengths = new ArrayList<Double>();
        final int[][] childBranch = new int[blockSize][];
        for (int i = 0; i < blockSize; i++) {
            NodeRef node = tree.getNode(postOrder[i]);
            children[i] = new int[tree.getChildCount(node)];
            childBranch[i] = new int[children[i].length];
            for (int j = 0; j < children[i].length; j++) {
                NodeRef child = tree.getChild(node, j);
                children[i][j] = child.getNumber();
                childBranch[i][j] = lengths.size();
                lengths.add(delegate.getBranchLength(child));
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
        final int topBranch = lengths.size();
        lengths.add(delegate.getBranchLength(tree.getNode(top)));
        final int branchCount = lengths.size();

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
        this.branchLengths = new double[branchCount];
        for (int e = 0; e < branchCount; e++) {
            branchLengths[e] = lengths.get(e);
        }
        this.logF = new double[branchCount][stateCount * stateCount];
        this.logFRate = Double.NaN;
        this.linearMatrices = new double[classWeights.length][][];
        ensureLogF(siteRateScale);

        final byte[][] blockStates = new byte[blockSize][];
        for (int i = 0; i < blockSize; i++) {
            blockStates[i] = sample ? store.getWritableStates(postOrder[i]) : store.getStates(postOrder[i]);
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

        for (int p = 0; p < patternCount; p++) {
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
                            blockStates[i][p] = -1;
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
                    blockStates[i][p] = (byte) (conditional.present[i] ? a : -1);
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
                        blockStates[i][p] = -1;
                    }
                }
                continue;
            }
            logProbability += choose(conditional.topLogP, conditional.topP, sample, blockStates[topIndex][p]);
            if (sample) {
                blockStates[topIndex][p] = (byte) chosenState;
            }
            // pre-order: each child given its (already chosen) parent
            for (int i = topIndex - 1; i >= 0; i--) {
                if (!conditional.present[i]) {
                    if (sample) {
                        blockStates[i][p] = -1;
                    }
                    continue;
                }
                int s = blockStates[parentIndex[i]][p];
                logProbability += choose(conditional.logRow(i, s), conditional.row(i, s), sample, blockStates[i][p]);
                if (sample) {
                    blockStates[i][p] = (byte) chosenState;
                }
            }
        }

        for (int n : postOrder) {
            inBlock[n] = false;
        }
        this.blockStore = null;
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
     */
    private final class Conditional {
        final double w;
        final double nu;
        final boolean[] present;
        final double[][] logPartials; // log density of everything below block node i given its state
        final int[] branchOf;         // the branch above each block node
        double[] topLogP;
        double[] topP;
        private final double[][][] logRows;
        private final double[][][] rows;

        Conditional(int p, double w, double nu) {
            final int blockSize = blockChildren.length;
            this.w = w;
            this.nu = nu;
            this.present = new boolean[blockSize];
            this.logPartials = new double[blockSize][stateCount];
            this.branchOf = new int[blockSize];
            this.logRows = new double[blockSize][stateCount][];
            this.rows = new double[blockSize][stateCount][];
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
                filterLog(p);
            }

            int topIndex = blockSize - 1;
            if (present[topIndex]) {
                double[] logWeights = new double[stateCount];
                if (blockTopParent == null) {
                    for (int s = 0; s < stateCount; s++) {
                        logWeights[s] = logPartials[topIndex][s] + w * logFrequencies[s];
                    }
                } else {
                    int parentState = blockStore.getState(blockTopParent.getNumber(), p);
                    if (parentState < 0) {
                        throw new IllegalStateException("Parent of the resampling block is ambiguous where the block is not");
                    }
                    double[] f = logF[blockTopBranch];
                    for (int s = 0; s < stateCount; s++) {
                        logWeights[s] = logPartials[topIndex][s] + w * f[parentState * stateCount + s];
                    }
                }
                topLogP = normaliseLog(logWeights);
                topP = exp(topLogP);
            }
        }

        /**
         * Filtering in linear space with exp(w logF) matrices (no exponentials per pattern).
         *
         * @return false if the partials underflow, in which case filterLog must be used
         */
        private boolean filterLinear(int p) {
            final double[][] m = getLinearMatrices(weightClass[p]);
            final int blockSize = blockChildren.length;
            final double[][] partials = new double[blockSize][stateCount];
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
                        int t = blockStore.getState(c, p);
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
                    for (int s = 0; s < stateCount; s++) {
                        partial[s] /= max;
                        logPartials[i][s] = Math.log(partial[s]);
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
            for (int i = 0; i < blockSize; i++) {
                double[] partial = logPartials[i];
                Arrays.fill(partial, 0.0);
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
                        int t = blockStore.getState(c, p);
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
         * @return the log distribution of (non-top) block node i given its parent has state s
         */
        double[] logRow(int i, int s) {
            if (logRows[i][s] == null) {
                ensureLogF(nu);
                double[] f = logF[branchOf[i]];
                double[] terms = new double[stateCount];
                for (int t = 0; t < stateCount; t++) {
                    terms[t] = w * f[s * stateCount + t] + logPartials[i][t];
                }
                logRows[i][s] = normaliseLog(terms);
                rows[i][s] = exp(logRows[i][s]);
            }
            return logRows[i][s];
        }

        double[] row(int i, int s) {
            logRow(i, s);
            return rows[i][s];
        }

        double[] logRowOrTop(int i, int s) {
            return i == present.length - 1 ? topLogP : logRow(i, s);
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

    private static double[] exp(double[] logP) {
        double[] p = new double[logP.length];
        for (int s = 0; s < logP.length; s++) {
            p[s] = Math.exp(logP[s]);
        }
        return p;
    }

    /**
     * For a pattern whose unambiguous nodes all have state a throughout the block's neighbourhood: the
     * probability that the block's unambiguous nodes all keep a, and for each node (in pre-order) the
     * probability that it is the first to differ from a.
     */
    private final class UniformConditional {
        final Conditional conditional;
        final int state;
        final double logProbabilityUnchanged;
        final double probabilityUnchanged;
        final double[] logPrefix;      // log probability that all nodes before this pre-order position keep a
        final double[] firstChangeMass; // probability that this position is the first to differ from a
        final double totalChangeMass;

        UniformConditional(Conditional conditional, int a) {
            this.conditional = conditional;
            this.state = a;
            int blockSize = conditional.present.length;
            logPrefix = new double[blockSize];
            firstChangeMass = new double[blockSize];
            double logPrefixSum = 0.0;
            double total = 0.0;
            for (int position = 0; position < blockSize; position++) {
                int i = blockSize - 1 - position;
                logPrefix[position] = logPrefixSum;
                if (!conditional.present[i]) {
                    // ambiguous nodes have no state to change
                    continue;
                }
                double[] logRow = conditional.logRowOrTop(i, a);
                double[] row = conditional.rowOrTop(i, a);
                double changed = 0.0;
                for (int t = 0; t < row.length; t++) {
                    if (t != a) {
                        changed += row[t];
                    }
                }
                firstChangeMass[position] = Math.exp(logPrefixSum) * changed;
                total += firstChangeMass[position];
                logPrefixSum += logRow[a];
            }
            logProbabilityUnchanged = logPrefixSum;
            probabilityUnchanged = Math.exp(logPrefixSum);
            totalChangeMass = total;
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

        double logProbability = uniform.logPrefix[first];

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
        blockStates[i][p] = (byte) newState;
        logProbability += logRow[newState];

        // later nodes in pre-order (whose parents are all already drawn) from their usual conditionals
        for (int position = first + 1; position < blockSize; position++) {
            int j = blockSize - 1 - position;
            if (!conditional.present[j]) {
                continue;
            }
            int parentState = blockStates[blockParentIndex[j]][p];
            logProbability += choose(conditional.logRow(j, parentState), conditional.row(j, parentState), true, -1);
            blockStates[j][p] = (byte) chosenState;
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
     * Log density of the path along a branch of length R (expected substitutions) from state a to
     * state b, allowing at most one substitution: -nu q_a R if a == b, otherwise
     * log(nu Q_ab R) - nu q_b R + log(expm1(x)/x) with x = nu (q_b - q_a) R.
     * Consistent with the likelihood computed by {@link MATDataLikelihoodDelegate}.
     */
    private void computeLogBranchDensities(double length, double nu, double[] out) {
        for (int a = 0; a < stateCount; a++) {
            for (int b = 0; b < stateCount; b++) {
                if (a == b) {
                    out[a * stateCount + b] = -nu * exitRates[a] * length;
                } else {
                    double rate = nu * infinitesimalMatrix[a * stateCount + b];
                    out[a * stateCount + b] = Math.log(rate * length) - nu * exitRates[b] * length
                            + MATDataLikelihoodDelegate.logExpm1OverX(nu * (exitRates[b] - exitRates[a]) * length);
                }
            }
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
