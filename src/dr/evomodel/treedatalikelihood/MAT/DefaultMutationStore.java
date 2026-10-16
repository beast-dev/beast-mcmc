/*
 * DefaultMutationStore.java
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

/**
 * A mutable mutation annotated tree with store/restore support.
 *
 * For each node it holds the states of all patterns and, for reconstruction by parsimony, the
 * Fitch state sets (as bitmasks, 0 meaning ambiguous). The mutations (including ambiguities) on the
 * branch above each node are derived from the states of the node and its parent by {@link #updateBranch}.
 *
 * The Fitch sets are only needed for (incremental) parsimony. Code that changes the topology and
 * states by other means marks them stale ({@link #setFitchSetsStale}) rather than keeping them up to
 * date, and parsimony recomputes them all before use.
 *
 * Each node's states carry a version number that changes whenever they may have changed, so that
 * quantities derived from them can be cached.
 *
 * Changes are journaled per node: the first change to a node after {@link #commit()} keeps the
 * previous arrays so that {@link #revert()} can put them back.
 */
public final class DefaultMutationStore implements MutationStore {

    private final int nodeCount;
    private final int patternCount;

    private int[] parents;
    private int root;
    private final byte[][] states;
    private final long[][] fitchSets;
    private boolean fitchSetsStale = false;
    private final BranchMutations[] mutations;
    private final long[] stateVersions;
    private long lastStateVersion = 0;

    // journal
    private int[] savedParents;
    private int savedRoot;
    private final byte[][] savedStates;
    private final long[][] savedFitchSets;
    private boolean savedFitchSetsStale;
    private boolean fitchSetsStaleChanged = false;
    private final BranchMutations[] savedMutations;
    private final boolean[] statesDirty;
    private final boolean[] fitchSetsDirty;
    private final boolean[] branchDirty;
    private final int[] dirtyNodes;
    private int dirtyNodeCount = 0;
    private final boolean[] nodeInJournal;

    public DefaultMutationStore(int[] parents, int patternCount) {
        this.nodeCount = parents.length;
        this.patternCount = patternCount;
        this.parents = parents.clone();
        this.root = findRoot(parents);
        this.states = new byte[nodeCount][patternCount];
        this.fitchSets = new long[nodeCount][patternCount];
        this.mutations = new BranchMutations[nodeCount];
        java.util.Arrays.fill(mutations, BranchMutations.EMPTY);
        this.stateVersions = new long[nodeCount];

        this.savedStates = new byte[nodeCount][];
        this.savedFitchSets = new long[nodeCount][];
        this.savedMutations = new BranchMutations[nodeCount];
        this.statesDirty = new boolean[nodeCount];
        this.fitchSetsDirty = new boolean[nodeCount];
        this.branchDirty = new boolean[nodeCount];
        this.nodeInJournal = new boolean[nodeCount];
        this.dirtyNodes = new int[nodeCount];
    }

    private static int findRoot(int[] parents) {
        for (int n = 0; n < parents.length; n++) {
            if (parents[n] == -1) {
                return n;
            }
        }
        throw new IllegalArgumentException("No root in parent array");
    }

    @Override
    public int getNodeCount() {
        return nodeCount;
    }

    @Override
    public int getPatternCount() {
        return patternCount;
    }

    @Override
    public int getParent(int nodeNumber) {
        return parents[nodeNumber];
    }

    @Override
    public int getRoot() {
        return root;
    }

    /**
     * @return true if the store describes the given topology
     */
    public boolean hasParents(int[] otherParents) {
        for (int n = 0; n < nodeCount; n++) {
            if (parents[n] != otherParents[n]) {
                return false;
            }
        }
        return true;
    }

    public void setParents(int[] newParents) {
        if (savedParents == null) {
            savedParents = parents;
            savedRoot = root;
        }
        parents = newParents.clone();
        root = findRoot(parents);
    }

    @Override
    public int getState(int nodeNumber, int pattern) {
        return states[nodeNumber][pattern];
    }

    /**
     * @return the node's states; callers must not modify the array (use {@link #getWritableStates})
     */
    public byte[] getStates(int nodeNumber) {
        return states[nodeNumber];
    }

    /**
     * @return the node's states for writing (each call counts as a change of the states)
     */
    public byte[] getWritableStates(int nodeNumber) {
        if (!statesDirty[nodeNumber]) {
            journal(nodeNumber);
            savedStates[nodeNumber] = states[nodeNumber];
            states[nodeNumber] = states[nodeNumber].clone();
            statesDirty[nodeNumber] = true;
        }
        stateVersions[nodeNumber] = ++lastStateVersion;
        return states[nodeNumber];
    }

    /**
     * @return a number that changes whenever the node's states may have changed
     */
    public long getStateVersion(int nodeNumber) {
        return stateVersions[nodeNumber];
    }

    public boolean isFitchSetsStale() {
        return fitchSetsStale;
    }

    public void setFitchSetsStale(boolean stale) {
        if (!fitchSetsStaleChanged) {
            savedFitchSetsStale = fitchSetsStale;
            fitchSetsStaleChanged = true;
        }
        fitchSetsStale = stale;
    }

    public long[] getFitchSets(int nodeNumber) {
        return fitchSets[nodeNumber];
    }

    public long[] getWritableFitchSets(int nodeNumber) {
        if (!fitchSetsDirty[nodeNumber]) {
            journal(nodeNumber);
            savedFitchSets[nodeNumber] = fitchSets[nodeNumber];
            fitchSets[nodeNumber] = fitchSets[nodeNumber].clone();
            fitchSetsDirty[nodeNumber] = true;
        }
        return fitchSets[nodeNumber];
    }

    @Override
    public BranchMutations getBranchMutations(int nodeNumber) {
        return mutations[nodeNumber];
    }

    @Override
    public int getMutationCount() {
        int count = 0;
        for (BranchMutations m : mutations) {
            count += m.getSubstitutionCount();
        }
        return count;
    }

    @Override
    public int getAmbiguityCount() {
        int count = 0;
        for (BranchMutations m : mutations) {
            count += m.getAmbiguityCount();
        }
        return count;
    }

    /**
     * Recomputes the mutations (including ambiguities) on the branch above the node from the states
     * of the node and its current parent.
     */
    public void updateBranch(int nodeNumber) {
        if (!branchDirty[nodeNumber]) {
            journal(nodeNumber);
            savedMutations[nodeNumber] = mutations[nodeNumber];
            branchDirty[nodeNumber] = true;
        }

        int parent = parents[nodeNumber];
        if (parent < 0) {
            mutations[nodeNumber] = BranchMutations.EMPTY;
            return;
        }

        byte[] childStates = states[nodeNumber];
        byte[] parentStates = states[parent];
        int substitutionCount = 0;
        int ambiguityCount = 0;
        for (int p = 0; p < patternCount; p++) {
            if (parentStates[p] >= 0 && childStates[p] != parentStates[p]) {
                if (childStates[p] < 0) {
                    ambiguityCount++;
                } else {
                    substitutionCount++;
                }
            }
        }
        if (substitutionCount + ambiguityCount == 0) {
            mutations[nodeNumber] = BranchMutations.EMPTY;
            return;
        }
        int[] patterns = new int[substitutionCount + ambiguityCount];
        byte[] from = new byte[patterns.length];
        byte[] to = new byte[patterns.length];
        int substitution = 0;
        int ambiguity = substitutionCount;
        for (int p = 0; p < patternCount; p++) {
            if (parentStates[p] >= 0 && childStates[p] != parentStates[p]) {
                int i = childStates[p] < 0 ? ambiguity++ : substitution++;
                patterns[i] = p;
                from[i] = parentStates[p];
                to[i] = childStates[p];
            }
        }
        mutations[nodeNumber] = new BranchMutations(patterns, from, to, substitutionCount);
    }

    /**
     * Updates the events on the branch above the node when its states, or its parent's, have changed only
     * at the given sites (and the parent is unchanged): events at other sites are kept.
     *
     * @param changed flags the sites in `sites`
     */
    public void updateBranchSites(int nodeNumber, boolean[] changed, int[] sites) {
        if (sites.length == 0) {
            return;
        }
        if (!branchDirty[nodeNumber]) {
            journal(nodeNumber);
            savedMutations[nodeNumber] = mutations[nodeNumber];
            branchDirty[nodeNumber] = true;
        }
        int parent = parents[nodeNumber];
        if (parent < 0) {
            mutations[nodeNumber] = BranchMutations.EMPTY;
            return;
        }
        BranchMutations old = mutations[nodeNumber];
        byte[] childStates = states[nodeNumber];
        byte[] parentStates = states[parent];

        int capacity = old.size() + sites.length;
        int[] subPatterns = new int[capacity];
        byte[] subFrom = new byte[capacity];
        byte[] subTo = new byte[capacity];
        int[] ambPatterns = new int[capacity];
        byte[] ambFrom = new byte[capacity];
        int subCount = 0;
        int ambCount = 0;
        for (int k = 0; k < old.size(); k++) {
            int p = old.getPattern(k);
            if (!changed[p]) {
                if (old.isAmbiguity(k)) {
                    ambPatterns[ambCount] = p;
                    ambFrom[ambCount++] = (byte) old.getFromState(k);
                } else {
                    subPatterns[subCount] = p;
                    subFrom[subCount] = (byte) old.getFromState(k);
                    subTo[subCount++] = (byte) old.getToState(k);
                }
            }
        }
        for (int p : sites) {
            if (parentStates[p] >= 0 && childStates[p] != parentStates[p]) {
                if (childStates[p] < 0) {
                    ambPatterns[ambCount] = p;
                    ambFrom[ambCount++] = parentStates[p];
                } else {
                    subPatterns[subCount] = p;
                    subFrom[subCount] = parentStates[p];
                    subTo[subCount++] = childStates[p];
                }
            }
        }
        if (subCount + ambCount == 0) {
            mutations[nodeNumber] = BranchMutations.EMPTY;
            return;
        }
        int[] patterns = new int[subCount + ambCount];
        byte[] from = new byte[patterns.length];
        byte[] to = new byte[patterns.length];
        System.arraycopy(subPatterns, 0, patterns, 0, subCount);
        System.arraycopy(subFrom, 0, from, 0, subCount);
        System.arraycopy(subTo, 0, to, 0, subCount);
        System.arraycopy(ambPatterns, 0, patterns, subCount, ambCount);
        System.arraycopy(ambFrom, 0, from, subCount, ambCount);
        java.util.Arrays.fill(to, subCount, patterns.length, (byte) -1);
        mutations[nodeNumber] = new BranchMutations(patterns, from, to, subCount);
    }

    /**
     * Replaces the events on the branch above the node (used by branch history models that sample
     * explicit histories).
     */
    public void setBranchMutations(int nodeNumber, BranchMutations branchMutations) {
        if (!branchDirty[nodeNumber]) {
            journal(nodeNumber);
            savedMutations[nodeNumber] = mutations[nodeNumber];
            branchDirty[nodeNumber] = true;
        }
        mutations[nodeNumber] = branchMutations;
    }

    public void updateAllBranches() {
        for (int n = 0; n < nodeCount; n++) {
            updateBranch(n);
        }
    }

    private void journal(int nodeNumber) {
        if (!nodeInJournal[nodeNumber]) {
            nodeInJournal[nodeNumber] = true;
            dirtyNodes[dirtyNodeCount++] = nodeNumber;
        }
    }

    /**
     * Accepts all changes since the last commit.
     */
    public void commit() {
        for (int i = 0; i < dirtyNodeCount; i++) {
            int n = dirtyNodes[i];
            savedStates[n] = null;
            savedFitchSets[n] = null;
            savedMutations[n] = null;
            statesDirty[n] = false;
            fitchSetsDirty[n] = false;
            branchDirty[n] = false;
            nodeInJournal[n] = false;
        }
        dirtyNodeCount = 0;
        savedParents = null;
        fitchSetsStaleChanged = false;
    }

    /**
     * Discards all changes since the last commit.
     */
    public void revert() {
        for (int i = 0; i < dirtyNodeCount; i++) {
            int n = dirtyNodes[i];
            if (statesDirty[n]) {
                states[n] = savedStates[n];
                stateVersions[n] = ++lastStateVersion;
            }
            if (fitchSetsDirty[n]) {
                fitchSets[n] = savedFitchSets[n];
            }
            if (branchDirty[n]) {
                mutations[n] = savedMutations[n];
            }
        }
        if (savedParents != null) {
            parents = savedParents;
            root = savedRoot;
        }
        if (fitchSetsStaleChanged) {
            fitchSetsStale = savedFitchSetsStale;
        }
        commit();
    }
}
