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

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Arrays;

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
 *
 * With discrete site rate categories ({@link #setCategories}), each variable pattern also has a category,
 * which is part of the sampled state, and each node a category composition: the summed weights of its
 * variable patterns in each (category, state), kept up to date as states and categories change. Both are
 * journaled.
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
    // site rate categories: the category of each pattern (-1 for constant patterns), or null without them
    private int categoryCount = 1;
    private byte[] categories;
    private double[][] categoryCompositions; // [node][category * stateCount + state]; null where not computed
    private int[] scratchPatterns;         // work arrays of updateBranch
    private byte[] scratchFrom;
    private byte[] scratchTo;

    // the nodes whose states or branch events have changed since the log was last cleared
    private final IndexSet changedNodes;

    // The journal: the first change to a node after commit() keeps what it replaces, so that revert() can
    // put it back. Per node (each with a flag saying whether it is in the journal): the parent, the states
    // with their version and composition, the Fitch sets and the branch events; then the root, the
    // categories and the stale flag.
    private final IndexSet journaledNodes;   // the nodes with anything in the journal
    private final boolean[] parentDirty;
    private final int[] savedParents;
    private final boolean[] statesDirty;
    private final byte[][] savedStates;
    private final long[] savedStateVersions;
    private final double[][] savedCompositions;
    private final boolean[] fitchSetsDirty;
    private final long[][] savedFitchSets;
    private final boolean[] branchDirty;
    private final BranchMutations[] savedMutations;
    private boolean rootSaved = false;
    private int savedRoot;
    private byte[] savedCategories;          // null unless a category has changed since the last commit
    private double[][] savedCategoryCompositions;
    private boolean[] categoryCompositionDirty;
    private boolean fitchSetsStaleChanged = false;
    private boolean savedFitchSetsStale;

    /**
     * @param patternWeights the weight of each pattern (for the compositions; may be longer than patternCount)
     * @param stateCount     the number of states
     */
    public DefaultMutationStore(int[] parents, int patternCount, double[] patternWeights, int stateCount) {
        this.nodeCount = parents.length;
        this.patternCount = patternCount;
        this.stateCount = stateCount;
        this.patternWeights = Arrays.copyOf(patternWeights, patternCount);
        this.parents = parents.clone();
        this.root = findRoot(parents);
        this.states = new byte[nodeCount][patternCount];
        this.fitchSets = new long[nodeCount][]; // allocated when parsimony first sets them
        this.mutations = new BranchMutations[nodeCount];
        Arrays.fill(mutations, BranchMutations.EMPTY);
        this.stateVersions = new long[nodeCount];
        this.compositions = new double[nodeCount][];
        this.changedNodes = new IndexSet(nodeCount);

        this.savedStates = new byte[nodeCount][];
        this.savedStateVersions = new long[nodeCount];
        this.savedCompositions = new double[nodeCount][];
        this.savedFitchSets = new long[nodeCount][];
        this.savedMutations = new BranchMutations[nodeCount];
        this.statesDirty = new boolean[nodeCount];
        this.fitchSetsDirty = new boolean[nodeCount];
        this.branchDirty = new boolean[nodeCount];
        this.journaledNodes = new IndexSet(nodeCount);
        this.savedParents = new int[nodeCount];
        this.parentDirty = new boolean[nodeCount];
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
            journaledNodes.add(nodeNumber);
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
        if (categories != null) {
            journalCategoryComposition(nodeNumber);
            categoryCompositions[nodeNumber] = null;
        }
        changedNodes.add(nodeNumber);
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
        if (categories != null && categories[pattern] >= 0 && categoryCompositions[nodeNumber] != null) {
            journalCategoryComposition(nodeNumber);
            double[] categoryComposition = categoryCompositions[nodeNumber];
            int offset = categories[pattern] * stateCount;
            double weight = patternWeights[pattern];
            if (old >= 0) {
                categoryComposition[offset + old] -= weight;
            }
            if (state >= 0) {
                categoryComposition[offset + state] += weight;
            }
        }
        stateVersions[nodeNumber] = ++lastStateVersion;
        changedNodes.add(nodeNumber);
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

    /**
     * Gives each pattern a site rate category (-1 for constant patterns, which have none). Not journaled:
     * for setting up the categories, before sampling.
     */
    public void setCategories(int categoryCount, byte[] patternCategories) {
        if (patternCategories.length != patternCount) {
            throw new IllegalArgumentException("Expected a category for each of the " + patternCount + " patterns");
        }
        this.categoryCount = categoryCount;
        this.categories = patternCategories.clone();
        this.categoryCompositions = new double[nodeCount][];
        this.savedCategoryCompositions = new double[nodeCount][];
        this.categoryCompositionDirty = new boolean[nodeCount];
        this.savedCategories = null;
    }

    /**
     * @return the number of site rate categories (1 without them)
     */
    public int getCategoryCount() {
        return categoryCount;
    }

    /**
     * @return the site rate category of the pattern (-1 for a constant pattern; 0 without categories)
     */
    public int getCategory(int pattern) {
        return categories == null ? 0 : categories[pattern];
    }

    /**
     * @return the categories of the patterns, or null without categories; callers must not modify it
     */
    public byte[] getCategories() {
        return categories;
    }

    /**
     * Changes the site rate category of a variable pattern, keeping the nodes' category compositions up
     * to date (in O(nodes)).
     */
    public void setCategory(int pattern, int category) {
        MATProfile.count(MATProfile.Counter.CATEGORY_CHANGES, categories[pattern] == category ? 0 : 1);
        final long profileStart = MATProfile.start();
        try {
            setCategoryImpl(pattern, category);
        } finally {
            MATProfile.stop(MATProfile.Phase.STORE_SET_CATEGORY, profileStart);
        }
    }

    private void setCategoryImpl(int pattern, int category) {
        final int old = categories[pattern];
        if (old == category) {
            return;
        }
        if (old < 0 || category < 0) {
            throw new IllegalArgumentException("Constant patterns have no site rate category");
        }
        if (savedCategories == null) {
            savedCategories = categories.clone();
        }
        categories[pattern] = (byte) category;
        final double weight = patternWeights[pattern];
        for (int n = 0; n < nodeCount; n++) {
            if (categoryCompositions[n] != null) {
                int state = states[n][pattern];
                if (state >= 0) {
                    journalCategoryComposition(n);
                    double[] categoryComposition = categoryCompositions[n];
                    categoryComposition[old * stateCount + state] -= weight;
                    categoryComposition[category * stateCount + state] += weight;
                }
            }
        }
    }

    /**
     * Replaces the site rate categories of all patterns (journaled), with the nodes' category compositions
     * for them, which the caller has computed (null entries are computed when next asked for). For changing
     * many categories at once: {@link #setCategory} updates every node for each pattern.
     */
    public void replaceCategories(byte[] patternCategories, double[][] nodeCategoryCompositions) {
        int changes = 0;
        for (int p = 0; p < patternCount; p++) {
            if (patternCategories[p] != categories[p]) {
                changes++;
            }
        }
        MATProfile.count(MATProfile.Counter.CATEGORY_CHANGES, changes);
        if (savedCategories == null) {
            savedCategories = categories;
        }
        categories = patternCategories.clone();
        for (int n = 0; n < nodeCount; n++) {
            journalCategoryComposition(n);
            categoryCompositions[n] = nodeCategoryCompositions[n];
        }
    }

    /**
     * @return the summed weights of the node's variable patterns in each category and state, indexed by
     * category * stateCount + state; callers must not modify it
     */
    public double[] getCategoryComposition(int nodeNumber) {
        double[] categoryComposition = categoryCompositions[nodeNumber];
        assert categoryComposition == null || Arrays.equals(categoryComposition, countCategoryComposition(nodeNumber)) :
                "category composition of node " + nodeNumber;
        if (categoryComposition == null) {
            if (statesDirty[nodeNumber] || savedCategories != null) {
                // computed from states or categories changed since the commit: revert must discard it
                journalCategoryComposition(nodeNumber);
            }
            categoryComposition = new double[categoryCount * stateCount];
            byte[] nodeStates = states[nodeNumber];
            for (int p = 0; p < patternCount; p++) {
                if (nodeStates[p] >= 0 && categories[p] >= 0) {
                    categoryComposition[categories[p] * stateCount + nodeStates[p]] += patternWeights[p];
                }
            }
            categoryCompositions[nodeNumber] = categoryComposition;
        }
        return categoryComposition;
    }

    private double[] countCategoryComposition(int nodeNumber) {
        double[] counts = new double[categoryCount * stateCount];
        byte[] nodeStates = states[nodeNumber];
        for (int p = 0; p < patternCount; p++) {
            if (nodeStates[p] >= 0 && categories[p] >= 0) {
                counts[categories[p] * stateCount + nodeStates[p]] += patternWeights[p];
            }
        }
        return counts;
    }

    private void journalCategoryComposition(int nodeNumber) {
        if (!categoryCompositionDirty[nodeNumber]) {
            journaledNodes.add(nodeNumber);
            savedCategoryCompositions[nodeNumber] = categoryCompositions[nodeNumber];
            if (categoryCompositions[nodeNumber] != null) {
                categoryCompositions[nodeNumber] = categoryCompositions[nodeNumber].clone();
            }
            categoryCompositionDirty[nodeNumber] = true;
        }
    }

    private void journalStates(int nodeNumber) {
        if (!statesDirty[nodeNumber]) {
            final long profileStart = MATProfile.start();
            journaledNodes.add(nodeNumber);
            savedStates[nodeNumber] = states[nodeNumber];
            states[nodeNumber] = states[nodeNumber].clone();
            savedStateVersions[nodeNumber] = stateVersions[nodeNumber];
            savedCompositions[nodeNumber] = compositions[nodeNumber];
            if (compositions[nodeNumber] != null) {
                compositions[nodeNumber] = compositions[nodeNumber].clone();
            }
            statesDirty[nodeNumber] = true;
            MATProfile.stop(MATProfile.Phase.STORE_JOURNAL_CLONE, profileStart);
            MATProfile.count(MATProfile.Counter.JOURNAL_CLONE_BYTES, patternCount);
        }
    }

    private void journalBranch(int nodeNumber) {
        if (!branchDirty[nodeNumber]) {
            journaledNodes.add(nodeNumber);
            savedMutations[nodeNumber] = mutations[nodeNumber];
            branchDirty[nodeNumber] = true;
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
            final long profileStart = MATProfile.start();
            journaledNodes.add(nodeNumber);
            long[] sets = fitchSets[nodeNumber];
            savedFitchSets[nodeNumber] = sets;
            fitchSets[nodeNumber] = sets == null ? new long[patternCount] : sets.clone();
            fitchSetsDirty[nodeNumber] = true;
            fitchSetsReleased = false;
            MATProfile.stop(MATProfile.Phase.STORE_FITCH_CLONE, profileStart);
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
        MATProfile.count(MATProfile.Counter.BRANCH_SCAN_PATTERNS, parents[nodeNumber] < 0 ? 0 : patternCount);
        final long profileStart = MATProfile.start();
        try {
            updateBranchImpl(nodeNumber);
        } finally {
            MATProfile.stop(MATProfile.Phase.STORE_UPDATE_BRANCH, profileStart);
        }
    }

    private void updateBranchImpl(int nodeNumber) {
        journalBranch(nodeNumber);
        changedNodes.add(nodeNumber);

        int parent = parents[nodeNumber];
        if (parent < 0) {
            mutations[nodeNumber] = BranchMutations.EMPTY;
            return;
        }

        final byte[] childStates = states[nodeNumber];
        final byte[] parentStates = states[parent];
        // views reading eight states at a time (to skip patterns where the two nodes' states are equal)
        final ByteBuffer childWords = ByteBuffer.wrap(childStates).order(ByteOrder.nativeOrder());
        final ByteBuffer parentWords = ByteBuffer.wrap(parentStates).order(ByteOrder.nativeOrder());
        ensureScratch();
        // one pass: substitutions from the front of the scratch arrays, ambiguities from the back (in reverse);
        // most patterns have the same state at both ends, which is skipped eight at a time
        int substitutionCount = 0;
        int ambiguityEnd = patternCount;
        int p = 0;
        final int wordEnd = patternCount - 7;
        while (p < patternCount) {
            if (p < wordEnd && childWords.getLong(p) == parentWords.getLong(p)) {
                p += 8;
                continue;
            }
            final int end = Math.min(p + 8, patternCount);
            for (; p < end; p++) {
                byte a = parentStates[p];
                byte b = childStates[p];
                if (a >= 0 && b != a) {
                    if (b < 0) {
                        ambiguityEnd--;
                        scratchPatterns[ambiguityEnd] = p;
                        scratchFrom[ambiguityEnd] = a;
                    } else {
                        scratchPatterns[substitutionCount] = p;
                        scratchFrom[substitutionCount] = a;
                        scratchTo[substitutionCount] = b;
                        substitutionCount++;
                    }
                }
            }
        }
        mutations[nodeNumber] = buildFromScratch(substitutionCount, ambiguityEnd);
    }

    private void ensureScratch() {
        if (scratchPatterns == null) {
            scratchPatterns = new int[patternCount];
            scratchFrom = new byte[patternCount];
            scratchTo = new byte[patternCount];
        }
    }

    /**
     * @return the events in the work arrays: substitutions in [0, substitutionCount), and ambiguities in
     * [ambiguityEnd, patternCount) in reverse pattern order
     */
    private BranchMutations buildFromScratch(int substitutionCount, int ambiguityEnd) {
        final int ambiguityCount = patternCount - ambiguityEnd;
        final int size = substitutionCount + ambiguityCount;
        if (size == 0) {
            return BranchMutations.EMPTY;
        }
        int[] patterns = new int[size];
        byte[] from = new byte[size];
        byte[] to = new byte[size];
        System.arraycopy(scratchPatterns, 0, patterns, 0, substitutionCount);
        System.arraycopy(scratchFrom, 0, from, 0, substitutionCount);
        System.arraycopy(scratchTo, 0, to, 0, substitutionCount);
        for (int i = 0; i < ambiguityCount; i++) {
            patterns[substitutionCount + i] = scratchPatterns[patternCount - 1 - i];
            from[substitutionCount + i] = scratchFrom[patternCount - 1 - i];
            to[substitutionCount + i] = -1;
        }
        return new BranchMutations(patterns, from, to, substitutionCount);
    }

    /**
     * Updates the events on the branch above the node when its states, or its parent's, have changed only
     * at the given sites (and the parent is unchanged): events at other sites are kept.
     *
     * @param changed the set of the sites in `sites`
     * @param sites   the sites, in pattern order (the events of every branch are kept in pattern order)
     */
    public void updateBranchSites(int nodeNumber, IndexSet changed, int[] sites) {
        MATProfile.count(MATProfile.Counter.BRANCH_SITE_UPDATES, sites.length);
        final long profileStart = MATProfile.start();
        try {
            updateBranchSitesImpl(nodeNumber, changed, sites);
        } finally {
            MATProfile.stop(MATProfile.Phase.STORE_UPDATE_BRANCH_SITES, profileStart);
        }
    }

    private void updateBranchSitesImpl(int nodeNumber, IndexSet changed, int[] sites) {
        if (sites.length == 0) {
            return;
        }
        journalBranch(nodeNumber);
        changedNodes.add(nodeNumber);
        int parent = parents[nodeNumber];
        if (parent < 0) {
            mutations[nodeNumber] = BranchMutations.EMPTY;
            return;
        }
        BranchMutations old = mutations[nodeNumber];
        byte[] childStates = states[nodeNumber];
        byte[] parentStates = states[parent];

        // The kept events (in pattern order) are merged with the new ones at the sites (in pattern order),
        // so that the events stay in pattern order. There is at most one event per pattern, so they fit in
        // the work arrays: substitutions from the front, ambiguities from the back (in reverse).
        ensureScratch();
        final int[] subPatterns = scratchPatterns;
        final byte[] subFrom = scratchFrom;
        final byte[] subTo = scratchTo;
        final int oldSubstitutionCount = old.getSubstitutionCount();
        final int oldSize = old.size();
        int subCount = 0;
        int ambEnd = patternCount;

        // substitutions: before each kept one, the new ones at the sites below its pattern
        int j = 0;
        for (int k = 0; k <= oldSubstitutionCount; k++) {
            int next = k < oldSubstitutionCount ? old.getPattern(k) : Integer.MAX_VALUE;
            for (; j < sites.length && sites[j] < next; j++) {
                int p = sites[j];
                if (parentStates[p] >= 0 && childStates[p] >= 0 && childStates[p] != parentStates[p]) {
                    subPatterns[subCount] = p;
                    subFrom[subCount] = parentStates[p];
                    subTo[subCount++] = childStates[p];
                }
            }
            if (k < oldSubstitutionCount && !changed.contains(next)) {
                subPatterns[subCount] = next;
                subFrom[subCount] = (byte) old.getFromState(k);
                subTo[subCount++] = (byte) old.getToState(k);
            }
        }

        // ambiguities likewise, from the back of the work arrays (so in reverse pattern order)
        j = 0;
        for (int k = oldSubstitutionCount; k <= oldSize; k++) {
            int next = k < oldSize ? old.getPattern(k) : Integer.MAX_VALUE;
            for (; j < sites.length && sites[j] < next; j++) {
                int p = sites[j];
                if (parentStates[p] >= 0 && childStates[p] < 0) {
                    ambEnd--;
                    subPatterns[ambEnd] = p;
                    subFrom[ambEnd] = parentStates[p];
                }
            }
            if (k < oldSize && !changed.contains(next)) {
                ambEnd--;
                subPatterns[ambEnd] = next;
                subFrom[ambEnd] = (byte) old.getFromState(k);
            }
        }
        mutations[nodeNumber] = buildFromScratch(subCount, ambEnd);
    }

    /**
     * Replaces the events on the branch above the node (used by branch history models that sample
     * explicit histories).
     */
    public void setBranchMutations(int nodeNumber, BranchMutations branchMutations) {
        journalBranch(nodeNumber);
        changedNodes.add(nodeNumber);
        mutations[nodeNumber] = branchMutations;
    }

    public void updateAllBranches() {
        for (int n = 0; n < nodeCount; n++) {
            updateBranch(n);
        }
    }

    /**
     * @return the number of nodes whose states or branch events have changed (including by revert)
     * since {@link #clearChangeLog}
     */
    public int getChangedNodeCount() {
        return changedNodes.size();
    }

    public int getChangedNode(int i) {
        return changedNodes.get(i);
    }

    public void clearChangeLog() {
        changedNodes.clear();
    }

    /**
     * Accepts all changes since the last commit.
     */
    public void commit() {
        final long profileStart = MATProfile.start();
        try {
            commitImpl();
        } finally {
            MATProfile.stop(MATProfile.Phase.STORE_COMMIT, profileStart);
        }
    }

    private void commitImpl() {
        for (int i = 0; i < journaledNodes.size(); i++) {
            int n = journaledNodes.get(i);
            savedStates[n] = null;
            savedCompositions[n] = null;
            savedFitchSets[n] = null;
            savedMutations[n] = null;
            statesDirty[n] = false;
            fitchSetsDirty[n] = false;
            branchDirty[n] = false;
            parentDirty[n] = false;
            if (categories != null) {
                savedCategoryCompositions[n] = null;
                categoryCompositionDirty[n] = false;
            }
        }
        journaledNodes.clear();
        rootSaved = false;
        savedCategories = null;
        fitchSetsStaleChanged = false;
        if (fitchSetsStale && !fitchSetsReleased) {
            // parsimony will recompute them all before use
            Arrays.fill(fitchSets, null);
            fitchSetsReleased = true;
        }
    }

    /**
     * Discards all changes since the last commit.
     */
    public void revert() {
        final long profileStart = MATProfile.start();
        try {
            revertImpl();
        } finally {
            MATProfile.stop(MATProfile.Phase.STORE_REVERT, profileStart);
        }
    }

    private void revertImpl() {
        for (int i = 0; i < journaledNodes.size(); i++) {
            int n = journaledNodes.get(i);
            if (statesDirty[n] || branchDirty[n]) {
                changedNodes.add(n);
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
            if (categories != null && categoryCompositionDirty[n]) {
                categoryCompositions[n] = savedCategoryCompositions[n];
            }
        }
        if (savedCategories != null) {
            categories = savedCategories;
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
