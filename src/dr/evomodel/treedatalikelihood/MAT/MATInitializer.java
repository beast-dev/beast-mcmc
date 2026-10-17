/*
 * MATInitializer.java
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
import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

/**
 * Builds and maintains a mutation annotated tree by Fitch parsimony.
 *
 * Any tip state that is not a single definite state (gaps, N, ?, and for simplicity also partial
 * ambiguities such as R) is treated as completely ambiguous: a pattern that is ambiguous in every tip
 * below a node gets state -1 there, which appears as an ambiguity on the branch where it becomes
 * ambiguous.
 *
 * The reconstruction is a deterministic function of the rooted topology, so the likelihood of the
 * resulting MAT is a well-defined function of the tree. {@link #update} recomputes it incrementally
 * after a topology change and gives exactly the same result as a full rebuild (provided the states
 * have not been changed by anything else, e.g. a {@link MutationHistoryResampler}).
 * For binary trees the result is a most-parsimonious reconstruction; multifurcations are handled
 * but are not guaranteed to be optimal.
 */
public class MATInitializer {

    private final PatternList patternList;
    private final Tree tree;
    private final int stateCount;
    private final int patternCount;
    private final long allStates;

    /**
     * @param tree the tree to reconstruct on; its topology may change between calls but its tips
     *             must keep their node numbers
     */
    public MATInitializer(PatternList patternList, Tree tree) {
        this.patternList = patternList;
        this.tree = tree;
        this.stateCount = patternList.getDataType().getStateCount();
        this.patternCount = patternList.getPatternCount();
        if (stateCount > 64) {
            throw new IllegalArgumentException("MAT likelihood supports at most 64 states, data type has " + stateCount);
        }
        this.allStates = stateCount == 64 ? -1L : (1L << stateCount) - 1;
    }

    public DefaultMutationStore buildMutationStore() {
        final int[] parents = getParents(tree);
        final int[][] children = getChildren(tree);
        final int[] postOrder = getPostOrder(tree);

        DefaultMutationStore store = new DefaultMutationStore(parents, patternCount, patternList.getPatternWeights(), stateCount);
        for (int n : postOrder) {
            if (children[n].length == 0) {
                System.arraycopy(getTipStateSets(tree.getNode(n)), 0, store.getWritableFitchSets(n), 0, patternCount);
            } else {
                updateFitchSets(store, n, children[n]);
            }
        }
        for (int i = postOrder.length - 1; i >= 0; i--) {
            updateStates(store, postOrder[i], parents[postOrder[i]]);
        }
        store.updateAllBranches();
        store.commit();
        return store;
    }

    /**
     * Brings the store up to date with the current topology of the tree, recomputing Fitch sets only
     * on paths from changed nodes to the root and states only where their inputs changed.
     *
     * @return the number of nodes whose states changed
     */
    public int update(DefaultMutationStore store) {
        final int nodeCount = tree.getNodeCount();
        // stale Fitch sets (the topology was changed without parsimony) force a full reconstruction
        final boolean full = store.isFitchSetsStale();
        if (full) {
            store.setFitchSetsStale(false);
        }
        final int[] parents = getParents(tree);
        final int[][] children = getChildren(tree);
        final int[] postOrder = getPostOrder(tree);
        if (full) {
            // the tips' sets follow from their states (partial ambiguities count as missing)
            for (int n = 0; n < nodeCount; n++) {
                if (children[n].length == 0) {
                    long[] sets = store.getWritableFitchSets(n);
                    for (int p = 0; p < patternCount; p++) {
                        int state = store.getState(n, p);
                        sets[p] = state < 0 ? 0 : 1L << state;
                    }
                }
            }
        }

        final boolean[] parentChanged = new boolean[nodeCount];
        final boolean[] needsSets = new boolean[nodeCount];
        for (int n = 0; n < nodeCount; n++) {
            int oldParent = store.getParent(n);
            if (parents[n] != oldParent) {
                parentChanged[n] = true;
                if (oldParent >= 0) {
                    needsSets[oldParent] = true;
                }
                if (parents[n] >= 0) {
                    needsSets[parents[n]] = true;
                }
            }
        }
        store.setParents(parents);

        final boolean[] setsChanged = new boolean[nodeCount];
        for (int n : postOrder) {
            if ((needsSets[n] || full) && children[n].length > 0 && (updateFitchSets(store, n, children[n]) || full)) {
                setsChanged[n] = true;
                if (parents[n] >= 0) {
                    needsSets[parents[n]] = true;
                }
            }
        }

        final boolean[] statesChanged = new boolean[nodeCount];
        int changedCount = 0;
        for (int i = postOrder.length - 1; i >= 0; i--) {
            int n = postOrder[i];
            int parent = parents[n];
            boolean parentStatesChanged = parent >= 0 && statesChanged[parent];
            if (setsChanged[n] || parentChanged[n] || parentStatesChanged) {
                statesChanged[n] = updateStates(store, n, parent);
                if (statesChanged[n]) {
                    changedCount++;
                }
            }
            if (parentChanged[n] || parentStatesChanged || statesChanged[n]) {
                store.updateBranch(n);
            }
        }
        return changedCount;
    }

    /**
     * Fitch upward step: a set of 0 means the pattern is ambiguous throughout the subtree.
     *
     * @return true if any set changed
     */
    boolean updateFitchSets(DefaultMutationStore store, int node, int[] children) {
        long[] current = store.getFitchSets(node);
        long[] sets = current == null ? store.getWritableFitchSets(node) : null;
        for (int p = 0; p < patternCount; p++) {
            long intersection = allStates;
            long union = 0;
            for (int c : children) {
                long childSet = store.getFitchSets(c)[p];
                if (childSet != 0) {
                    intersection &= childSet;
                    union |= childSet;
                }
            }
            long set = union == 0 ? 0 : (intersection != 0 ? intersection : union);
            if (current == null || set != current[p]) {
                if (sets == null) {
                    sets = store.getWritableFitchSets(node);
                }
                sets[p] = set;
            }
        }
        return sets != null;
    }

    /**
     * Fitch downward step: keep the parental state wherever the node's set allows it.
     *
     * @return true if any state changed
     */
    private boolean updateStates(DefaultMutationStore store, int node, int parent) {
        long[] sets = store.getFitchSets(node);
        byte[] current = store.getStates(node);
        byte[] parentStates = parent >= 0 ? store.getStates(parent) : null;
        boolean changed = false;
        for (int p = 0; p < patternCount; p++) {
            long set = sets[p];
            int state;
            if (set == 0) {
                state = -1;
            } else if (parentStates != null && parentStates[p] >= 0 && ((set >>> parentStates[p]) & 1L) != 0) {
                state = parentStates[p];
            } else {
                state = Long.numberOfTrailingZeros(set);
            }
            if (state != current[p]) {
                store.setState(node, p, state);
                changed = true;
            }
        }
        return changed;
    }

    /**
     * @return the tip's state sets as bitmasks by pattern, with any ambiguous state as 0
     */
    private long[] getTipStateSets(NodeRef tip) {
        DataType dataType = patternList.getDataType();
        String taxonId = tree.getNodeTaxon(tip).getId();
        int taxonIndex = patternList.getTaxonIndex(taxonId);
        if (taxonIndex < 0) {
            throw new IllegalArgumentException("Taxon " + taxonId + " is not in the pattern list");
        }
        long[] sets = new long[patternCount];
        for (int p = 0; p < patternCount; p++) {
            boolean[] stateSet = dataType.getStateSet(patternList.getPatternState(taxonIndex, p));
            long set = 0;
            for (int s = 0; s < stateCount; s++) {
                if (stateSet[s]) {
                    set |= 1L << s;
                }
            }
            // partial ambiguities are treated as complete ambiguities
            sets[p] = Long.bitCount(set) == 1 ? set : 0;
        }
        return sets;
    }

    /**
     * @return node numbers in post-order (children before parents), computed without recursion
     */
    static int[] getPostOrder(Tree tree) {
        int nodeCount = tree.getNodeCount();
        int[] order = new int[nodeCount];
        int[] stack = new int[nodeCount];
        int top = 0;
        int k = nodeCount;
        stack[top++] = tree.getRoot().getNumber();
        // fill from the end in pre-order (parent before children) so that read forwards it is a post-order
        while (top > 0) {
            int n = stack[--top];
            order[--k] = n;
            NodeRef node = tree.getNode(n);
            for (int i = 0; i < tree.getChildCount(node); i++) {
                stack[top++] = tree.getChild(node, i).getNumber();
            }
        }
        return order;
    }

    /**
     * @return the child node numbers of each node
     */
    static int[][] getChildren(Tree tree) {
        int[][] children = new int[tree.getNodeCount()][];
        for (int n = 0; n < children.length; n++) {
            NodeRef node = tree.getNode(n);
            children[n] = new int[tree.getChildCount(node)];
            for (int i = 0; i < children[n].length; i++) {
                children[n][i] = tree.getChild(node, i).getNumber();
            }
        }
        return children;
    }

    /**
     * @return the parent node number of each node, -1 for the root
     */
    static int[] getParents(Tree tree) {
        int[] parents = new int[tree.getNodeCount()];
        for (int n = 0; n < parents.length; n++) {
            NodeRef parent = tree.getParent(tree.getNode(n));
            parents[n] = parent == null ? -1 : parent.getNumber();
        }
        return parents;
    }

    public PatternList getPatternList() {
        return patternList;
    }

    public Tree getTree() {
        return tree;
    }
}
