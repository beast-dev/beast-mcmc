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
 * date, and parsimony recomputes them all before use. Stale sets are released when committed (they are
 * eight times the size of the states), and allocated again when parsimony next needs them.
 *
 * Each node's states carry a version number that changes whenever they may have changed, so that
 * quantities derived from them can be cached. Each node also has a summary of its states, its
 * composition (the summed pattern weights of its sites in each state), which {@link #setState} keeps up
 * to date as states are changed.
 *
 * Changes are journaled per node: the first change to a node after {@link #commit()} keeps the
 * previous arrays (and the version and composition of the states) so that {@link #revert()} can put
 * them back.
 */
public final class DefaultMutationStore implements MutationStore {

    private final int nodeCount;
    private final int patternCount;
    private final int stateCount;
    private final double[] patternWeights;

    private final int[] parents;
    private int root;
    private final byte[][] states;
    private final long[][] fitchSets;
    private boolean fitchSetsStale = false;
    private boolean fitchSetsReleased = false;
    private final BranchMutations[] mutations;
    private final long[] stateVersions;
    private long lastStateVersion = 0;
    private final double[][] compositions; // null where not yet computed from the states

    // the nodes whose states or branch events have changed since the log was last cleared
    private final boolean[] changeLogged;
    private final int[] changeLog;
    private int changeLogCount = 0;

    // journal
    private final int[] savedParents;
    private final boolean[] parentDirty;
    private int savedRoot;
    private boolean rootSaved = false;
    private final byte[][] savedStates;
    private final long[] savedStateVersions;
    private final double[][] savedCompositions;
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

    /**
     * @param patternWeights the weight of each pattern (for the compositions; may be longer than patternCount)
     * @param stateCount     the number of states
     */
    public DefaultMutationStore(int[] parents, int patternCount, double[] patternWeights, int stateCount) {
        this.nodeCount = parents.length;
        this.patternCount = patternCount;
        this.stateCount = stateCount;
        this.patternWeights = java.util.Arrays.copyOf(patternWeights, patternCount);
        this.parents = parents.clone();
        this.root = findRoot(parents);
        this.states = new byte[nodeCount][patternCount];
        this.fitchSets = new long[nodeCount][]; // allocated when parsimony first sets them
        this.mutations = new BranchMutations[nodeCount];
        java.util.Arrays.fill(mutations, BranchMutations.EMPTY);
        this.stateVersions = new long[nodeCount];
        this.compositions = new double[nodeCount][];
        this.changeLogged = new boolean[nodeCount];
        this.changeLog = new int[nodeCount];

        this.savedStates = new byte[nodeCount][];
        this.savedStateVersions = new long[nodeCount];
        this.savedCompositions = new double[nodeCount][];
        this.savedFitchSets = new long[nodeCount][];
        this.savedMutations = new BranchMutations[nodeCount];
        this.statesDirty = new boolean[nodeCount];
        this.fitchSetsDirty = new boolean[nodeCount];
        this.branchDirty = new boolean[nodeCount];
        this.nodeInJournal = new boolean[nodeCount];
        this.savedParents = new int[nodeCount];
        this.parentDirty = new boolean[nodeCount];
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
        for (int n = 0; n < nodeCount; n++) {
            setParent(n, newParents[n]);
        }
        if (root < 0 || parents[root] >= 0) {
            root = findRoot(parents);
        }
    }

    /**
     * Sets the parent of a node (-1 makes it the root). A topology change must leave exactly one node
     * without a parent.
     */
    public void setParent(int nodeNumber, int parent) {
        if (parents[nodeNumber] == parent) {
            return;
        }
        if (!rootSaved) {
            savedRoot = root;
            rootSaved = true;
        }
        if (!parentDirty[nodeNumber]) {
            journal(nodeNumber);
            savedParents[nodeNumber] = parents[nodeNumber];
            parentDirty[nodeNumber] = true;
        }
        parents[nodeNumber] = parent;
        if (parent < 0) {
            root = nodeNumber;
        } else if (root == nodeNumber) {
            root = -1; // until the new root's parent is cleared
        }
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
     * @return the node's states for writing directly (each call counts as a change of the states, and the
     * node's composition is recomputed in full when next asked for); to change a few states, use
     * {@link #setState}
     */
    public byte[] getWritableStates(int nodeNumber) {
        journalStates(nodeNumber);
        stateVersions[nodeNumber] = ++lastStateVersion;
        compositions[nodeNumber] = null;
        logChange(nodeNumber);
        return states[nodeNumber];
    }

    /**
     * @return the node's states, as the array that {@link #setState} will change until the next commit or
     * revert; callers must not modify it themselves
     */
    public byte[] getStatesForUpdate(int nodeNumber) {
        journalStates(nodeNumber);
        return states[nodeNumber];
    }

    /**
     * Sets the state of a node at a pattern (-1 for ambiguous), keeping its composition up to date.
     */
    public void setState(int nodeNumber, int pattern, int state) {
        int old = states[nodeNumber][pattern];
        if (old == state) {
            return;
        }
        journalStates(nodeNumber);
        states[nodeNumber][pattern] = (byte) state;
        double[] composition = compositions[nodeNumber];
        if (composition != null) {
            double weight = patternWeights[pattern];
            if (old >= 0) {
                composition[old] -= weight;
            }
            if (state >= 0) {
                composition[state] += weight;
            }
        }
        stateVersions[nodeNumber] = ++lastStateVersion;
        logChange(nodeNumber);
    }

    /**
     * @return the summed pattern weights of the node's sites in each state; callers must not modify it
     */
    public double[] getComposition(int nodeNumber) {
        double[] composition = compositions[nodeNumber];
        if (composition == null) {
            composition = new double[stateCount];
            byte[] nodeStates = states[nodeNumber];
            for (int p = 0; p < patternCount; p++) {
                if (nodeStates[p] >= 0) {
                    composition[nodeStates[p]] += patternWeights[p];
                }
            }
            compositions[nodeNumber] = composition;
        }
        return composition;
    }

    private void journalStates(int nodeNumber) {
        if (!statesDirty[nodeNumber]) {
            journal(nodeNumber);
            savedStates[nodeNumber] = states[nodeNumber];
            states[nodeNumber] = states[nodeNumber].clone();
            savedStateVersions[nodeNumber] = stateVersions[nodeNumber];
            savedCompositions[nodeNumber] = compositions[nodeNumber];
            if (compositions[nodeNumber] != null) {
                compositions[nodeNumber] = compositions[nodeNumber].clone();
            }
            statesDirty[nodeNumber] = true;
        }
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

    /**
     * @return the node's Fitch sets, or null if they have not been computed (or were released when stale)
     */
    public long[] getFitchSets(int nodeNumber) {
        return fitchSets[nodeNumber];
    }

    public long[] getWritableFitchSets(int nodeNumber) {
        if (!fitchSetsDirty[nodeNumber]) {
            journal(nodeNumber);
            long[] sets = fitchSets[nodeNumber];
            savedFitchSets[nodeNumber] = sets;
            fitchSets[nodeNumber] = sets == null ? new long[patternCount] : sets.clone();
            fitchSetsDirty[nodeNumber] = true;
            fitchSetsReleased = false;
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
        logChange(nodeNumber);

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
     * @param sites   the sites, in pattern order (the events of every branch are kept in pattern order)
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
        logChange(nodeNumber);
        int parent = parents[nodeNumber];
        if (parent < 0) {
            mutations[nodeNumber] = BranchMutations.EMPTY;
            return;
        }
        BranchMutations old = mutations[nodeNumber];
        byte[] childStates = states[nodeNumber];
        byte[] parentStates = states[parent];

        // the kept events (in pattern order) are merged with the new ones at the sites (in pattern order),
        // so that the events stay in pattern order
        int capacity = old.size() + sites.length;
        int[] subPatterns = new int[capacity];
        byte[] subFrom = new byte[capacity];
        byte[] subTo = new byte[capacity];
        int[] ambPatterns = new int[capacity];
        byte[] ambFrom = new byte[capacity];
        int subCount = 0;
        int ambCount = 0;
        int j = 0;
        for (int k = 0; k <= old.getSubstitutionCount(); k++) {
            int next = k < old.getSubstitutionCount() ? old.getPattern(k) : Integer.MAX_VALUE;
            for (; j < sites.length && sites[j] < next; j++) {
                int p = sites[j];
                if (parentStates[p] >= 0 && childStates[p] >= 0 && childStates[p] != parentStates[p]) {
                    subPatterns[subCount] = p;
                    subFrom[subCount] = parentStates[p];
                    subTo[subCount++] = childStates[p];
                }
            }
            if (k < old.getSubstitutionCount() && !changed[next]) {
                subPatterns[subCount] = next;
                subFrom[subCount] = (byte) old.getFromState(k);
                subTo[subCount++] = (byte) old.getToState(k);
            }
        }
        j = 0;
        for (int k = old.getSubstitutionCount(); k <= old.size(); k++) {
            int next = k < old.size() ? old.getPattern(k) : Integer.MAX_VALUE;
            for (; j < sites.length && sites[j] < next; j++) {
                int p = sites[j];
                if (parentStates[p] >= 0 && childStates[p] < 0) {
                    ambPatterns[ambCount] = p;
                    ambFrom[ambCount++] = parentStates[p];
                }
            }
            if (k < old.size() && !changed[next]) {
                ambPatterns[ambCount] = next;
                ambFrom[ambCount++] = (byte) old.getFromState(k);
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
        logChange(nodeNumber);
        mutations[nodeNumber] = branchMutations;
    }

    public void updateAllBranches() {
        for (int n = 0; n < nodeCount; n++) {
            updateBranch(n);
        }
    }

    private void logChange(int nodeNumber) {
        if (!changeLogged[nodeNumber]) {
            changeLogged[nodeNumber] = true;
            changeLog[changeLogCount++] = nodeNumber;
        }
    }

    /**
     * @return the number of nodes whose states or branch events have changed (including by revert)
     * since {@link #clearChangeLog}
     */
    public int getChangedNodeCount() {
        return changeLogCount;
    }

    public int getChangedNode(int i) {
        return changeLog[i];
    }

    public void clearChangeLog() {
        for (int i = 0; i < changeLogCount; i++) {
            changeLogged[changeLog[i]] = false;
        }
        changeLogCount = 0;
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
            savedCompositions[n] = null;
            savedFitchSets[n] = null;
            savedMutations[n] = null;
            statesDirty[n] = false;
            fitchSetsDirty[n] = false;
            branchDirty[n] = false;
            parentDirty[n] = false;
            nodeInJournal[n] = false;
        }
        dirtyNodeCount = 0;
        rootSaved = false;
        fitchSetsStaleChanged = false;
        if (fitchSetsStale && !fitchSetsReleased) {
            // parsimony will recompute them all before use
            java.util.Arrays.fill(fitchSets, null);
            fitchSetsReleased = true;
        }
    }

    /**
     * Discards all changes since the last commit.
     */
    public void revert() {
        for (int i = 0; i < dirtyNodeCount; i++) {
            int n = dirtyNodes[i];
            if (statesDirty[n] || branchDirty[n]) {
                logChange(n);
            }
            if (statesDirty[n]) {
                // the states are as they were at the saved version, so caches keyed by it are valid again
                states[n] = savedStates[n];
                stateVersions[n] = savedStateVersions[n];
                compositions[n] = savedCompositions[n];
            }
            if (fitchSetsDirty[n]) {
                fitchSets[n] = savedFitchSets[n];
            }
            if (branchDirty[n]) {
                mutations[n] = savedMutations[n];
            }
            if (parentDirty[n]) {
                parents[n] = savedParents[n];
            }
        }
        if (rootSaved) {
            root = savedRoot;
        }
        if (fitchSetsStaleChanged) {
            fitchSetsStale = savedFitchSetsStale;
        }
        commit();
    }
}
