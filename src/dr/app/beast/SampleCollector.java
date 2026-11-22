/*
 * SampleCollector.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
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
 *
 */

package dr.app.beast;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeLogger;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Logger;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.NumberColumn;
import dr.inference.markovchain.MarkovChain;
import dr.inference.markovchain.MarkovChainListener;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCCriterion;
import dr.inference.model.Model;

import java.util.*;

/**
 * Keeps the samples of one chain in memory for convergence diagnostics: the chain's posterior,
 * the numeric columns of its main log, and the frequencies of the clades of its logged tree.
 *
 * Trees are not stored. Each clade is identified by a 64-bit hash (the sum of random keys of
 * its taxa, which needs O(n) per tree) and counted in blocks of consecutive samples, so a
 * burn-in can be dropped later with memory that does not grow with the chain length.
 *
 * For a parallel tempering set it listens to every member and records only the one at
 * temperature 1.
 *
 * @author Andrew Rambaut
 */
public class SampleCollector implements MarkovChainListener {

    public static final int TREE_BLOCK_COUNT = 50;

    /**
     * @param name  a label for reports
     * @param mcmcs the members of a chain set, with the first one built with the real loggers
     */
    public SampleCollector(String name, MCMC[] mcmcs) {
        this.name = name;

        // The main log is the one with the most numeric columns (screen logs have fewer)
        // and the first tree logger gives the tree.
        int[] bestLogger = new int[mcmcs.length];
        int[] treeLogger = new int[mcmcs.length];
        Arrays.fill(bestLogger, -1);
        Arrays.fill(treeLogger, -1);

        Logger[] loggers0 = mcmcs[0].getLoggers();
        int mostColumns = 0;
        for (int j = 0; loggers0 != null && j < loggers0.length; j++) {
            if (loggers0[j] instanceof TreeLogger) {
                if (treeLogger[0] < 0) {
                    treeLogger[0] = j;
                }
            } else if (loggers0[j] instanceof MCLogger) {
                int count = getNumberColumns((MCLogger) loggers0[j]).size();
                if (count > mostColumns) {
                    mostColumns = count;
                    bestLogger[0] = j;
                }
            }
        }

        // every member was built from the same XML so the loggers are in the same order
        valueColumns = new List[mcmcs.length];
        trees = new Tree[mcmcs.length];
        String[] labels = new String[0];
        long logEvery = 0;
        long treeEvery = 0;
        for (int i = 0; i < mcmcs.length; i++) {
            Logger[] loggers = mcmcs[i].getLoggers();
            if (bestLogger[0] >= 0) {
                MCLogger logger = (MCLogger) loggers[bestLogger[0]];
                valueColumns[i] = getNumberColumns(logger);
                logEvery = logger.getLogEvery();
            } else {
                valueColumns[i] = Collections.emptyList();
            }
            if (treeLogger[0] >= 0) {
                TreeLogger logger = (TreeLogger) loggers[treeLogger[0]];
                trees[i] = logger.getTree();
                treeEvery = logger.getLogEvery();
            }
        }

        labels = new String[valueColumns[0].size()];
        for (int j = 0; j < labels.length; j++) {
            labels[j] = valueColumns[0].get(j).getLabel();
        }
        this.labels = labels;
        this.treeEvery = (trees[0] != null ? treeEvery : 0);

        long chainLength = mcmcs[0].getChainLength();
        // the posterior is always sampled, at the main log's interval if there is one
        if (logEvery <= 0) {
            logEvery = (this.treeEvery > 0 ? this.treeEvery : Math.max(1, chainLength / 1000));
        }
        this.logEvery = logEvery;
        long expectedTreeSamples = (this.treeEvery > 0 ? chainLength / this.treeEvery : 0);
        this.treeSamplesPerBlock = Math.max(1, expectedTreeSamples / TREE_BLOCK_COUNT);

        this.chainIndex = new IdentityHashMap<MarkovChain, Integer>();
        for (int i = 0; i < mcmcs.length; i++) {
            chainIndex.put(mcmcs[i].getMarkovChain(), i);
        }
        this.coldChain = mcmcs[0].getMarkovChain();
    }

    private static List<NumberColumn> getNumberColumns(MCLogger logger) {
        List<NumberColumn> columns = new ArrayList<NumberColumn>();
        for (int i = 0; i < logger.getColumnCount(); i++) {
            LogColumn column = logger.getColumn(i);
            if (column instanceof NumberColumn) {
                columns.add((NumberColumn) column);
            }
        }
        return columns;
    }

    /**
     * Starts listening to the chains. Must be called before the chains start.
     */
    public void attach(MCMC[] mcmcs) {
        for (MCMC mcmc : mcmcs) {
            mcmc.getMarkovChain().addMarkovChainListener(this);
        }
    }

    public String getName() {
        return name;
    }

    public String[] getLabels() {
        return labels;
    }

    public boolean hasTrees() {
        return treeEvery > 0;
    }

    // MarkovChainListener: called on the chain's own thread

    @Override
    public void currentState(long state, MarkovChain markovChain, Model currentModel) {
        boolean logValues = logEvery > 0 && state % logEvery == 0;
        boolean logTree = treeEvery > 0 && state % treeEvery == 0;
        if (!logValues && !logTree) {
            return;
        }
        if (((MCMCCriterion) markovChain.getAcceptor()).getTemperature() != 1.0) {
            return;
        }
        int i = chainIndex.get(markovChain);
        coldChain = markovChain;

        double score = markovChain.getCurrentScore();
        double[] row = null;
        if (logValues) {
            List<NumberColumn> columns = valueColumns[i];
            row = new double[columns.size()];
            for (int j = 0; j < row.length; j++) {
                row[j] = columns.get(j).getDoubleValue();
            }
        }
        long[] clades = (logTree ? getCladeKeys(trees[i]) : null);

        synchronized (this) {
            lastState = state;
            if (row != null) {
                rows.add(row);
                scores.add(score);
            }
            if (clades != null) {
                if (currentBlock == null || currentBlockSize >= treeSamplesPerBlock) {
                    currentBlock = new HashMap<Long, Integer>();
                    currentBlockSize = 0;
                    treeBlocks.add(currentBlock);
                    treeBlockSizes.add(0);
                }
                for (long key : clades) {
                    Integer count = currentBlock.get(key);
                    currentBlock.put(key, count == null ? 1 : count + 1);
                }
                currentBlockSize += 1;
                treeBlockSizes.set(treeBlockSizes.size() - 1, currentBlockSize);
            }
        }
    }

    @Override
    public void bestState(long state, MarkovChain markovChain, Model bestModel) {
    }

    @Override
    public void finished(long chainLength, MarkovChain markovChain) {
    }

    /**
     * @return the keys of the non-trivial clades of the tree (excluding the tips and the root)
     */
    private long[] getCladeKeys(Tree tree) {
        int nodeCount = tree.getNodeCount();
        int tipCount = tree.getExternalNodeCount();
        long[] keys = new long[nodeCount];
        long[] clades = new long[Math.max(0, tipCount - 2)];
        int cladeCount = 0;

        // iterative post-order traversal (the trees can be deep)
        NodeRef root = tree.getRoot();
        Deque<NodeRef> stack = new ArrayDeque<NodeRef>();
        Deque<NodeRef> order = new ArrayDeque<NodeRef>();
        stack.push(root);
        while (!stack.isEmpty()) {
            NodeRef node = stack.pop();
            order.push(node);
            for (int j = 0; j < tree.getChildCount(node); j++) {
                stack.push(tree.getChild(node, j));
            }
        }
        while (!order.isEmpty()) {
            NodeRef node = order.pop();
            int n = node.getNumber();
            if (tree.isExternal(node)) {
                keys[n] = getTaxonKey(tree.getNodeTaxon(node).getId());
            } else {
                long key = 0;
                for (int j = 0; j < tree.getChildCount(node); j++) {
                    key += keys[tree.getChild(node, j).getNumber()];
                }
                keys[n] = key;
                if (node != root && cladeCount < clades.length) {
                    clades[cladeCount++] = key;
                }
            }
        }
        return (cladeCount == clades.length ? clades : Arrays.copyOf(clades, cladeCount));
    }

    private long getTaxonKey(String id) {
        Long key = taxonKeys.get(id);
        if (key == null) {
            // FNV-1a then a SplitMix64 finaliser: the same for the same id in every chain
            long h = 0xcbf29ce484222325L;
            for (int i = 0; i < id.length(); i++) {
                h ^= id.charAt(i);
                h *= 0x100000001b3L;
            }
            h = (h ^ (h >>> 30)) * 0xbf58476d1ce4e5b9L;
            h = (h ^ (h >>> 27)) * 0x94d049bb133111ebL;
            key = h ^ (h >>> 31);
            taxonKeys.put(id, key);
        }
        return key;
    }

    /**
     * A copy of the samples collected so far, safe to use on another thread.
     */
    public synchronized Snapshot snapshot() {
        List<Map<Long, Integer>> blocks = new ArrayList<Map<Long, Integer>>(treeBlocks);
        if (!blocks.isEmpty()) {
            // only the last block is still being added to
            int last = blocks.size() - 1;
            blocks.set(last, new HashMap<Long, Integer>(blocks.get(last)));
        }
        return new Snapshot(lastState, coldChain.getCurrentScore(), new ArrayList<double[]>(rows),
                new ArrayList<Double>(scores), blocks, new ArrayList<Integer>(treeBlockSizes), logEvery);
    }

    public static class Snapshot {
        Snapshot(long state, double posterior, List<double[]> rows, List<Double> scores,
                 List<Map<Long, Integer>> treeBlocks, List<Integer> treeBlockSizes, long logEvery) {
            this.state = state;
            this.posterior = posterior;
            this.rows = rows;
            this.scores = scores;
            this.treeBlocks = treeBlocks;
            this.treeBlockSizes = treeBlockSizes;
            this.logEvery = logEvery;
        }

        public final long state;
        public final double posterior;
        public final List<double[]> rows;
        /**
         * the chain's log posterior at the same states as rows
         */
        public final List<Double> scores;
        public final List<Map<Long, Integer>> treeBlocks;
        public final List<Integer> treeBlockSizes;
        public final long logEvery;
    }

    private final String name;
    private final String[] labels;
    private final List<NumberColumn>[] valueColumns;
    private final Tree[] trees;
    private final long logEvery;
    private final long treeEvery;
    private final long treeSamplesPerBlock;
    private final Map<MarkovChain, Integer> chainIndex;

    // only used on the chain threads (one at a time for a tempering set)
    private final Map<String, Long> taxonKeys = new HashMap<String, Long>();

    // guarded by this
    private final List<double[]> rows = new ArrayList<double[]>();
    private final List<Double> scores = new ArrayList<Double>();
    private final List<Map<Long, Integer>> treeBlocks = new ArrayList<Map<Long, Integer>>();
    private final List<Integer> treeBlockSizes = new ArrayList<Integer>();
    private Map<Long, Integer> currentBlock = null;
    private int currentBlockSize = 0;
    private long lastState = 0;
    private volatile MarkovChain coldChain;
}
