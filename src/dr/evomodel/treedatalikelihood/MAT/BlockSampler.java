/*
 * BlockSampler.java
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
import dr.math.MathUtils;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Scores or draws the states of a block of a MAT's internal nodes, one site at a time, from their full
 * conditional given all the other states in the current tree (the block part of
 * {@link MutationHistoryResampler}): the {@link Block} of one call, its {@link Conditional}s (forward
 * filtering, backward sampling) and the conditionals of the sites whose state is the same throughout the
 * block's neighbourhood ({@link UniformConditional}).
 */
final class BlockSampler {

    private final MATDataLikelihoodDelegate delegate;
    private final DefaultMutationStore store;
    private final BranchHistoryModel historyModel;
    private final Tree tree;
    private final int stateCount;
    private final double[] patternWeights;
    // patterns grouped by weight, so conditionals can be looked up by [state][weight class]
    private final int[] weightClass;
    private final double[] classWeights;
    private final int defaultWeightClass;   // of the variable patterns, when sites are handled sparsely
    // the sites' rate classes: their site rate categories (0 with a single category), or none (-1) with
    // per-site rates
    private final int rateClassCount;
    // the nodes of the current block (see Block), and each one's index in it
    private final boolean[] inBlock;
    private final int[] nodeBlockIndex;

    BlockSampler(MATDataLikelihoodDelegate delegate, int[] weightClass, double[] classWeights,
                 int defaultWeightClass, int rateClassCount) {
        this.delegate = delegate;
        this.store = delegate.getMutableMutationStore();
        this.historyModel = delegate.getHistoryModel();
        this.tree = delegate.getTree();
        this.stateCount = delegate.getPatternList().getDataType().getStateCount();
        this.patternWeights = delegate.getPatternWeights();
        this.weightClass = weightClass;
        this.classWeights = classWeights;
        this.defaultWeightClass = defaultWeightClass;
        this.rateClassCount = rateClassCount;
        this.inBlock = new boolean[tree.getNodeCount()];
        this.nodeBlockIndex = new int[tree.getNodeCount()];
    }

    /**
     * @return the block of the given nodes in the current tree, with the store's states of its nodes (for
     * writing if sampling); its nodes are marked in inBlock until it is released
     */
    Block newBlock(int[] blockNodes, boolean sample) {
        return new Block(blockNodes, sample);
    }

    /**
     * Scores or draws the block's states at a variable site, by the site's class around the block.
     *
     * @return logProbability plus the log probability of the existing (or drawn) states under their full
     * conditional (the terms are added to the running total one at a time, as the Hastings ratio sums them)
     */
    double sampleSite(Block block, int p, byte patternClass, double logProbability) {
        final double w = patternWeights[p];
        final double nu = delegate.getPatternRate(p);
        if (patternClass == MutationHistoryResampler.NO_EVENTS) {
            return logProbability + sampleUniformSite(block, p, w, nu);
        }
        final ConditionalKey key = block.keyFor(p, w, nu);
        final Conditional conditional = block.conditionalFor(key, p, w, nu);
        if (patternClass == MutationHistoryResampler.AMBIGUITIES_ONLY) {
            return logProbability + sampleAmbiguityOnlySite(block, p, key, conditional);
        }
        return sampleSubstitutionSite(block, p, conditional, logProbability);
    }

    /**
     * A site with no events around the block: every node in and around it has the same state.
     */
    private double sampleUniformSite(Block block, int p, double w, double nu) {
        final int a = block.states[0][p]; // (the states are not yet overwritten)
        if (a < 0) {
            return 0.0;
        }
        final UniformConditional uniform = block.uniformConditionalFor(p, a, w, nu);
        if (!block.sample || MathUtils.nextDouble() < uniform.probabilityUnchanged) {
            return uniform.logProbabilityUnchanged;
        }
        return sampleChangedUniformPattern(block, uniform, p);
    }

    /**
     * A site whose only events around the block are ambiguities: every unambiguous node in and around the
     * block has the same state, before and after the move.
     */
    private double sampleAmbiguityOnlySite(Block block, int p, ConditionalKey key, Conditional conditional) {
        final int blockSize = block.nodes.length;
        int a = -1;
        for (int n : block.neighbours) {
            a = Math.max(a, store.getState(n, p));
        }
        for (int i = 0; i < blockSize; i++) {
            a = Math.max(a, block.states[i][p]);
        }
        if (a < 0 || !conditional.present[block.topIndex]) {
            block.clearStates(p);
            return 0.0;
        }
        final UniformConditional uniform = block.uniformConditionalFor(key, conditional, a);
        if (!block.sample) {
            return uniform.logProbabilityUnchanged;
        }
        // the unchanged configuration: unambiguous nodes take a (which they may not have had before the move)
        for (int i = 0; i < blockSize; i++) {
            block.setState(i, p, (byte) (conditional.present[i] ? a : -1));
        }
        if (MathUtils.nextDouble() < uniform.probabilityUnchanged) {
            return uniform.logProbabilityUnchanged;
        }
        return sampleChangedUniformPattern(block, uniform, p);
    }

    /**
     * A site with substitutions around the block: the top node from its marginal, then each other node
     * given its (already chosen) parent, in pre-order.
     *
     * @return logProbability plus the log probability of each node's state, added in that order
     */
    private double sampleSubstitutionSite(Block block, int p, Conditional conditional, double logProbability) {
        final int topIndex = block.topIndex;
        if (!conditional.present[topIndex]) {
            block.clearStates(p);
            return logProbability;
        }
        final byte[][] blockStates = block.states;
        final boolean sample = block.sample;
        final int top = choose(conditional.topP(), sample, blockStates[topIndex][p]);
        logProbability += logProbabilityOf(conditional.topLogP(), top);
        if (sample) {
            block.setState(topIndex, p, (byte) top);
        }
        for (int i = topIndex - 1; i >= 0; i--) {
            if (!conditional.present[i]) {
                if (sample) {
                    block.setState(i, p, -1);
                }
                continue;
            }
            final int s = blockStates[block.parentIndex[i]][p];
            final int state = choose(conditional.row(i, s), sample, blockStates[i][p]);
            logProbability += logProbabilityOf(conditional.logRow(i, s), state);
            if (sample) {
                block.setState(i, p, (byte) state);
            }
        }
        return logProbability;
    }

    /**
     * One call of sampleBlock: the block's nodes and branches in the current tree, the log branch densities
     * (for the rate of the patterns being handled), the block's states, which are only changed through
     * {@link #setState} so that the store's node summaries stay up to date, and the conditionals computed
     * so far. The block's nodes are marked in inBlock until {@link #release}.
     */
    final class Block {
        final boolean sample;
        final int[] nodes;          // in post-order, the top last
        final int topIndex;
        final NodeRef topParent;
        final int[][] children;     // the children of each block node
        final int[][] childBranch;  // the index of the branch above each child
        final int[] index;          // node number -> block index (only read for block nodes)
        final int[] parentIndex;
        final int[] neighbours;     // the nodes outside the block whose branches touch it: the top node's
                                    // parent (if any), then the non-block children of each block node
        final int topBranch;        // the index of the branch above the top node
        final double[] branchLengths;
        final double[][] logF;      // the log branch densities, for the rate logFRate (with per-site rates)
        double logFRate = Double.NaN;
        final double[][][] classLogF;          // the log branch densities of each rate class
        final double[][][][] linearMatrices;   // [weight class][rate class]
        final double[] logFrequencies;
        final byte[][] states;

        // The full conditional of a pattern's block states depends only on its weight, its rate and the
        // states adjacent to the block, so it is computed once for each distinct combination of these.
        // Draws are still made independently for each pattern.
        private final Map<ConditionalKey, Conditional> cache = new HashMap<ConditionalKey, Conditional>();
        private final Map<ConditionalKey, UniformConditional> uniformCache =
                new HashMap<ConditionalKey, UniformConditional>();
        // without per-site rates, conditionals for patterns with no events depend only on (state, weight, rate class)
        private final UniformConditional[][] uniformByClass;
        private final byte[] boundary; // work array of keyFor

        Block(int[] blockNodes, boolean sample) {
            this.sample = sample;
            this.nodes = getBlockPostOrder(blockNodes);
            final int blockSize = nodes.length;
            this.topIndex = blockSize - 1;
            this.topParent = tree.getParent(tree.getNode(nodes[topIndex]));
            for (int n : nodes) {
                inBlock[n] = true;
            }

            this.children = new int[blockSize][];
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

            int neighbourCount = topParent == null ? 0 : 1;
            for (int i = 0; i < blockSize; i++) {
                for (int c : children[i]) {
                    if (!inBlock[c]) {
                        neighbourCount++;
                    }
                }
            }
            this.neighbours = new int[neighbourCount];
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
            this.boundary = new byte[topParent == null ? neighbourCount + 1 : neighbourCount];

            double[] frequencies = delegate.getRootFrequencyModel().getFrequencies();
            this.logFrequencies = new double[stateCount];
            for (int a = 0; a < stateCount; a++) {
                logFrequencies[a] = Math.log(frequencies[a]);
            }
            this.logF = new double[branchTotal][stateCount * stateCount];
            this.classLogF = new double[rateClassCount][][];
            this.linearMatrices = new double[classWeights.length][rateClassCount][][];
            this.uniformByClass = delegate.hasSiteRates() ? null :
                    new UniformConditional[stateCount][classWeights.length * rateClassCount];

            this.states = new byte[blockSize][];
            for (int i = 0; i < blockSize; i++) {
                states[i] = sample ? store.getStatesForUpdate(nodes[i]) : store.getStates(nodes[i]);
            }
        }

        /**
         * @return the key of the site's conditional: the states around the block (-2 for a missing parent),
         * the weight and the rate
         */
        ConditionalKey keyFor(int p, double w, double nu) {
            int k = 0;
            if (topParent == null) {
                boundary[k++] = -2;
            }
            for (int n : neighbours) {
                boundary[k++] = (byte) store.getState(n, p);
            }
            return new ConditionalKey(boundary.clone(), w, nu);
        }

        /**
         * @return the conditional of the site's block states, shared with the sites that have the same key
         */
        Conditional conditionalFor(ConditionalKey key, int p, double w, double nu) {
            Conditional conditional = cache.get(key);
            if (conditional == null) {
                conditional = siteConditional(this, p, w, nu);
                cache.put(key, conditional);
            }
            return conditional;
        }

        /**
         * @return the conditional of a site with the given key whose unambiguous nodes all have state a
         */
        UniformConditional uniformConditionalFor(ConditionalKey key, Conditional conditional, int a) {
            UniformConditional uniform = uniformCache.get(key);
            if (uniform == null) {
                uniform = new UniformConditional(conditional, a);
                uniformCache.put(key, uniform);
            }
            return uniform;
        }

        /**
         * @return the conditional of a site with no events around the block, with state a throughout
         */
        UniformConditional uniformConditionalFor(int p, int a, double w, double nu) {
            if (uniformByClass != null) {
                final int c = weightClass[p] * rateClassCount + rateClassOf(p);
                UniformConditional uniform = uniformByClass[a][c];
                if (uniform == null) {
                    uniform = new UniformConditional(siteConditional(this, p, w, nu), a);
                    uniformByClass[a][c] = uniform;
                }
                return uniform;
            }
            final ConditionalKey key = new ConditionalKey(new byte[]{(byte) a}, w, nu);
            UniformConditional uniform = uniformCache.get(key);
            if (uniform == null) {
                uniform = new UniformConditional(siteConditional(this, p, w, nu), a);
                uniformCache.put(key, uniform);
            }
            return uniform;
        }

        /**
         * @return the log branch densities for a rate class (or, for rate class -1, the given rate)
         */
        double[][] getLogF(int rateClass, double nu) {
            if (rateClass < 0) {
                if (nu != logFRate) {
                    for (int e = 0; e < logF.length; e++) {
                        historyModel.getLogTransitionDensities(branchLengths[e], nu, logF[e]);
                    }
                    logFRate = nu;
                }
                return logF;
            }
            double[][] f = classLogF[rateClass];
            if (f == null) {
                double rate = rateOfClass(rateClass);
                f = new double[logF.length][stateCount * stateCount];
                for (int e = 0; e < f.length; e++) {
                    historyModel.getLogTransitionDensities(branchLengths[e], rate, f[e]);
                }
                classLogF[rateClass] = f;
            }
            return f;
        }

        /**
         * @return exp(w logF) for every branch, for the weight class's w and the rate class
         */
        double[][] getLinearMatrices(int weightClassIndex, int rateClass) {
            double[][] matrices = linearMatrices[weightClassIndex][rateClass];
            if (matrices == null) {
                double w = classWeights[weightClassIndex];
                double[][] f = getLogF(rateClass, rateOfClass(rateClass));
                matrices = new double[f.length][stateCount * stateCount];
                for (int e = 0; e < f.length; e++) {
                    for (int j = 0; j < stateCount * stateCount; j++) {
                        matrices[e][j] = Math.exp(w * f[e][j]);
                    }
                }
                linearMatrices[weightClassIndex][rateClass] = matrices;
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

        /**
         * Makes the site ambiguous at every block node (when drawing; scoring leaves the states alone).
         */
        void clearStates(int p) {
            if (sample) {
                for (int i = 0; i < states.length; i++) {
                    setState(i, p, -1);
                }
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

    /**
     * @return the rate class of a pattern: its site rate category, or -1 with per-site rates
     */
    private int rateClassOf(int p) {
        return delegate.hasSiteRates() ? -1 : store.getCategory(p);
    }

    /**
     * @return the rate of a rate class
     */
    double rateOfClass(int rateClass) {
        return rateClassCount > 1 ? delegate.getCategoryRates()[rateClass] : delegate.getSiteRateScale();
    }


    /**
     * The full conditional distribution of a block's states for one pattern, factorised as the
     * distribution of the top node and of each other node given its parent. The upward (filtering)
     * pass is done on construction; the distribution of a node given its parent's state is computed
     * when first needed.
     * <p>
     * Without per-site rates the filtering is done in linear space and the distributions are
     * normalised in linear space; their logarithms are only taken when asked for.
     */
    final class Conditional {
        final Block block;
        final double w;
        final double nu;
        private final int rateClass;
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

        /**
         * @param p            the pattern whose states around the block the conditional is for, or -1
         * @param uniformState if p is -1, the state at every node around the block
         */
        private Conditional(Block block, int p, int uniformState, double w, double nu, int rateClass) {
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

            this.rateClass = rateClass;
            // (with per-site rates, rate class -1, there are too many rates for linear matrices)
            if (rateClass < 0 || !filterLinear(p, rateClass)) {
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
                        double[] f = logF()[block.topBranch];
                        for (int s = 0; s < stateCount; s++) {
                            logWeights[s] = logPartials[topIndex][s] + w * f[parentState * stateCount + s];
                        }
                    }
                    topLogP = normaliseLog(logWeights);
                    topP = exp(topLogP);
                }
            }
        }

        /**
         * @return the log branch densities at its rate (with per-site rates the block keeps those of one rate
         * at a time, so they are brought back to this one)
         */
        private double[][] logF() {
            return block.getLogF(rateClass, nu);
        }

        private int stateOf(int node, int p) {
            return uniformState >= 0 ? uniformState : store.getState(node, p);
        }

        /**
         * Filtering in linear space with exp(w logF) matrices (no exponentials or logarithms per pattern).
         *
         * @return false if the partials underflow, in which case filterLog must be used
         */
        private boolean filterLinear(int p, int rateClass) {
            final double[][] m = block.getLinearMatrices(p >= 0 ? weightClass[p] : defaultWeightClass, rateClass);
            final int[][] children = block.children;
            final int[][] childBranch = block.childBranch;
            final int[] index = block.index;
            final boolean[] inBlock = BlockSampler.this.inBlock;
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
            final double[][] logF = logF();
            final boolean[] inBlock = BlockSampler.this.inBlock;
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
                    double[] f = logF()[branchOf[i]];
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
     * @return the full conditional of the block's states at a pattern, given the states around the block
     */
    Conditional siteConditional(Block block, int p, double w, double nu) {
        return new Conditional(block, p, -1, w, nu, rateClassOf(p));
    }

    /**
     * @return the full conditional of the block's states at a site of the rate class that has state a at
     * every node around the block
     */
    Conditional uniformStateConditional(Block block, int a, double w, double nu, int rateClass) {
        return new Conditional(block, -1, a, w, nu, rateClass);
    }

    /**
     * For a pattern whose unambiguous nodes all have state a throughout the block's neighbourhood: the
     * probability that the block's unambiguous nodes all keep a, and for each node (in pre-order) the
     * probability that it is the first to differ from a. Computed in linear space with one logarithm.
     */
    static final class UniformConditional {
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
    double sampleChangedUniformPattern(Block block, UniformConditional uniform, int p) {
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
            final int parentState = block.states[block.parentIndex[j]][p];
            final int state = choose(conditional.row(j, parentState), true, -1);
            logProbability += logProbabilityOf(conditional.logRow(j, parentState), state);
            block.setState(j, p, (byte) state);
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
     * Draws a state from a normalised distribution or, if not sampling, takes the existing state.
     *
     * @return the state; 0 if the distribution has no mass (the proposal is then impossible, see
     * {@link #logProbabilityOf})
     */
    private static int choose(double[] p, boolean sample, int state) {
        if (!sample) {
            if (state < 0) {
                throw new IllegalStateException("MAT has an ambiguous state at a node that has data below it");
            }
            return state;
        }
        if (Double.isNaN(p[0])) {
            // the conditional has no mass: the proposal is impossible
            return 0;
        }
        double u = MathUtils.nextDouble();
        double cumulative = 0.0;
        for (int s = 0; s < p.length; s++) {
            cumulative += p[s];
            if (u < cumulative) {
                return s;
            }
        }
        // rounding: take the last state with non-zero probability
        for (int s = p.length - 1; s >= 0; s--) {
            if (p[s] > 0.0) {
                return s;
            }
        }
        return 0;
    }

    /**
     * @return the log probability of the state under a normalised log distribution: -infinity where the
     * distribution has no mass (NaN throughout, see normaliseLinear)
     */
    private static double logProbabilityOf(double[] logP, int state) {
        return Double.isNaN(logP[state]) ? Double.NEGATIVE_INFINITY : logP[state];
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
